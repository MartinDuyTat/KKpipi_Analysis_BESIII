// Martin Duy Tat 3rd October 2022

#include<vector>
#include<memory>
#include<stdexcept>
#include<utility>
#include"TMatrixTSym.h"
#include"TRandom.h"
#include"TFile.h"
#include"TMatrixTSym.h"
#include"TMath.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedDTYields::RawBinnedDTYields(const std::vector<AsymmetricUncertainty> &Yields,
				     const TMatrixTSym<double> &CorrelationMatrix):
  m_Yields(Yields),
  m_CorrelationMatrix(CorrelationMatrix) {
}

const std::vector<AsymmetricUncertainty>&
RawBinnedDTYields::GetDoubleTagYields() const {
  return m_Yields;
}

const TMatrixTSym<double>& RawBinnedDTYields::GetCorrelationMatrix() const {
  return m_CorrelationMatrix;
}

std::pair<std::vector<AsymmetricUncertainty>, double>
RawBinnedDTYields::GetToyYields(const std::vector<double> &PredictedYields,
				bool SymmetricUncertainties) const {
  if(m_Yields.size() != PredictedYields.size()) {
    throw std::runtime_error("Yields and predicted yields not the same size vector");
  }
  double Probability = 1.0;
  auto ToyYields = m_Yields;
  for(std::size_t i = 0; i < m_Yields.size(); i++) {
    const auto GeneratedYield = GenerateYield(PredictedYields[i],
					      m_Yields[i],
					      SymmetricUncertainties);
    ToyYields[i].Value = GeneratedYield.first;
    /*if(!SymmetricUncertainties) {
      const auto Uncertainties = GetAsymmetricUncertainties(GeneratedYield.first);
      ToyYields[i].PlusUncertainty = Uncertainties.first;
      ToyYields[i].NegativeUncertainty = Uncertainties.second;
    }*/
    Probability *= GeneratedYield.second;
  }
  return std::make_pair(ToyYields, Probability);
}

std::pair<double, double> RawBinnedDTYields::GenerateYield(
  double PredictedYield,
  const AsymmetricUncertainty &DataYield,
  bool SymmetricUncertainties) const {
  if(SymmetricUncertainties) {
    // Symmetric uncertainties are straight forward, use Gaussian
    const double Sigma = DataYield.SymmetricUncertainty;
    const double GeneratedYield = gRandom->Gaus(PredictedYield, Sigma);
    const double GeneratedProbability = TMath::Gaus(GeneratedYield,
						    PredictedYield,
						    Sigma);
    return std::make_pair(GeneratedYield, GeneratedProbability);
  } else {
    // Asymmetric uncertainties are more complicated
    const double SigmaPlus = DataYield.PlusUncertainty;
    const double SigmaMinus = DataYield.NegativeUncertainty;
    const double Std2 = SigmaPlus*SigmaMinus;
    const double StdDiff = SigmaMinus - SigmaPlus;
    const double SigmaFactor = SigmaPlus*SigmaMinus
                              /TMath::Abs(SigmaMinus - SigmaPlus);
    // There are some hard limits where the likelihood diverges
    double LowerLimit = PredictedYield;
    double UpperLimit = PredictedYield;
    if(SigmaPlus > SigmaMinus) {
      LowerLimit -= 10.0*SigmaPlus;
      UpperLimit += SigmaFactor;
    } else if(SigmaMinus > SigmaPlus) {
      LowerLimit -= SigmaFactor;
      UpperLimit += 10.0*SigmaMinus;
    } else {
      LowerLimit -= 5.0*SigmaMinus;
      UpperLimit += 5.0*SigmaPlus;
    }
    LowerLimit = std::max(LowerLimit, 0.0);
    double GeneratedProbability, GeneratedYield;
    // Perform a rejection sampling using the asymmetric likelihood
    while(true) {
      GeneratedYield = gRandom->Uniform(LowerLimit, UpperLimit);
      const double YieldDiff = GeneratedYield - PredictedYield;
      double Var = Std2 + StdDiff*YieldDiff;
      GeneratedProbability = TMath::Exp(-YieldDiff*YieldDiff/(2.0*Var));
      const double RandomNumber = gRandom->Uniform(0, 1);
      if(RandomNumber < GeneratedProbability) {
	break;
      }
    }
    return std::make_pair(GeneratedYield, GeneratedProbability);
  }
}

TMatrixTSym<double>
RawBinnedDTYields::LoadCorrelationMatrix(const std::string &Tag,
					    const Settings &settings) const {
  const std::string Filename =
    Utilities::ReplaceString(settings.get("RawYields_CorrMatrix"), "TAG", Tag);
  return LoadCorrelationMatrix(Filename);
}

TMatrixTSym<double>
RawBinnedDTYields::LoadCorrelationMatrix(const std::string &Filename) const {
  TFile File(Filename.c_str());
  TMatrixTSym<double> *CorrelationMatrix = nullptr;
  File.GetObject("CorrelationMatrix", CorrelationMatrix);
  File.Close();
  return *CorrelationMatrix;
}
