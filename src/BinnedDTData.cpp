// Martin Duy Tat 9th October 2022

#include<vector>
#include<memory>
#include<algorithm>
#include<stdexcept>
#include<iostream>
#include<iomanip>
#include<fstream>
#include"Eigen/Dense"
#include"TMatrixT.h"
#include"TMath.h"
#include"TRandom.h"
#include"TF1.h"
#include"Math/Functor.h"
#include"Math/RootFinder.h"
#include"BinnedDTData.h"
#include"RawBinnedDTYields.h"
#include"RawBinnedCPTagYields.h"
#include"RawBinnedSCMBTagYields.h"
#include"RawBinnedFlavourTagYields.h"
#include"BinnedDTYieldPrediction.h"
#include"BinnedCPTagYieldPrediction.h"
#include"BinnedSCMBTagYieldPrediction.h"
#include"BinnedFlavourTagYieldPrediction.h"
#include"RawBinnedDTYieldLikelihood.h"
#include"cisiFitterParameters.h"

BinnedDTData::BinnedDTData(const std::string &Tag,
			   const Settings &settings):
  m_TagMode(Tag),
  m_FullLikelihood(settings.contains(Tag + "_FullLikelihood") &&
		   settings.getB(Tag + "_FullLikelihood")),
  m_DTYields(GetRawDTYields(Tag, settings)),
  m_DTYieldLikelihood(GetDTYieldLikelihood(Tag, settings)),
  m_DTPredictions(GetDTPredictions(Tag, settings)),
  m_SymmetricUncertainties(settings.getB("SymmetricUncertainties")),
  m_Settings(settings) {
}

double BinnedDTData::GetLogLikelihood(
  const cisiFitterParameters &Parameters) const {
  if(m_FullLikelihood) {
    const auto PredictedYields =
      m_DTPredictions->GetPredictedBinYields(Parameters);
    return m_DTYieldLikelihood->GetLogLikelihood(PredictedYields);
  } else {
    const auto MeasuredYields = m_DTYields->GetDoubleTagYields();
    return GetLogLikelihood(Parameters, MeasuredYields);
  }
}

double BinnedDTData::GetLogLikelihood(
  const cisiFitterParameters &Parameters,
  const std::vector<AsymmetricUncertainty> &MeasuredYields) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(Parameters);
  const auto CorrelationMatrix = m_DTYields->GetCorrelationMatrix();
  const auto InvCovMatrix = CreateInvCovarianceMatrix(PredictedYields,
						      MeasuredYields,
						      CorrelationMatrix);
  double LogLikelihood = 0.0;
  std::size_t Size = MeasuredYields.size();
  for(std::size_t i = 0; i < Size; i++) {
    const double Diff_i = MeasuredYields[i].Value - PredictedYields[i];
    for(std::size_t j = 0; j <= i; j++) {
      const double Diff_j = MeasuredYields[j].Value - PredictedYields[j];
      const int SymmetryFactor = (i == j ? 1 : 2);
      LogLikelihood += InvCovMatrix(i, j)*Diff_i*Diff_j*SymmetryFactor;
    }
  }
  return LogLikelihood;
}

void BinnedDTData::LoadToyDataset(int ToyNumber) const {
  m_DTYields = GetRawDTYields(m_TagMode, m_Settings, ToyNumber);
  if(m_FullLikelihood) {
    m_DTYieldLikelihood = GetDTYieldLikelihood(m_TagMode, m_Settings, ToyNumber);
  }
}

void BinnedDTData::PrintComparison(const cisiFitterParameters &Parameters) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(Parameters);
  const auto MeasuredYields = m_DTYields->GetDoubleTagYields();
  const double LL = GetLogLikelihood(Parameters, MeasuredYields);
  std::cout << "Fitted and predicted yield comparison for " << m_TagMode << ":\n";
  std::cout << "Log-likelihood = " << LL << "\n";
  std::cout << std::left << std::setw(10) << "Fitted";
  std::cout << std::left << std::setw(10) << "+ Error";
  std::cout << std::left << std::setw(10) << "- Error";
  std::cout << std::left << std::setw(10) << "Predicted\n";
  for(std::size_t i = 0; i < PredictedYields.size(); i++) {
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << MeasuredYields[i].Value;
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << MeasuredYields[i].PlusUncertainty;
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << MeasuredYields[i].NegativeUncertainty;
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << PredictedYields[i] << "\n";
  }
}

double BinnedDTData::GetAsymmetricStd(const AsymmetricUncertainty &Yield,
				      double Prediction) const {
  if(m_SymmetricUncertainties) {
    return Yield.SymmetricUncertainty;
  }
  double NegativeUncertainty = Yield.NegativeUncertainty;
  if(NegativeUncertainty <= 0.0 && Yield.PlusUncertainty != 0.0) {
    // Just set it to something close to zero to avoid numerical issues
    NegativeUncertainty = 0.2;
  }
  const double Std2 = Yield.PlusUncertainty*NegativeUncertainty;
  const double StdDiff = NegativeUncertainty - Yield.PlusUncertainty;
  const double YieldDiff = Yield.Value - Prediction;
  const double Var = Std2 + StdDiff*YieldDiff;
  if(Var < 0.0) {
    return 1.0e-6;
  } else {
    return TMath::Sqrt(Var);
  }
}

Eigen::MatrixXd BinnedDTData::CreateInvCovarianceMatrix(
  const std::vector<double> &PredictedYields,
  const std::vector<AsymmetricUncertainty> &MeasuredYields,
  const TMatrixT<double> &CorrelationMatrix) const {
  std::size_t Size = PredictedYields.size();
  Eigen::MatrixXd CovMatrix(Size, Size);
  for(std::size_t i = 0; i < Size; i++) {
    const double Std_i = GetAsymmetricStd(MeasuredYields[i], PredictedYields[i]);
    for(std::size_t j = 0; j <= i; j++) {
      const double Std_j = GetAsymmetricStd(MeasuredYields[j], PredictedYields[j]);
      CovMatrix(i, j) = CorrelationMatrix(i, j)*Std_i*Std_j;
      if(i != j) {
	CovMatrix(j, i) = CovMatrix(i, j);
      }
    }
  }
  return CovMatrix.llt().solve(Eigen::MatrixXd::Identity(Size, Size));
}

std::unique_ptr<const RawBinnedDTYields> BinnedDTData::GetRawDTYields(
  const std::string &Tag,
  const Settings &settings,
  int ToyNumber) const {
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), Tag);
  const auto iter_Flavour = std::find(m_FlavourTags.begin(),
				      m_FlavourTags.end(), Tag);
  if(iter_CP != m_CPTags.end()) {
    return std::make_unique<const RawBinnedCPTagYields>(Tag,
							settings,
							ToyNumber);
  } else if(iter_SCMB != m_SCMBTags.end()) {
    return std::make_unique<const RawBinnedSCMBTagYields>(Tag,
							  settings,
							  ToyNumber);
  } else if(iter_Flavour != m_FlavourTags.end()) {
    return std::make_unique<const RawBinnedFlavourTagYields>(Tag,
							     settings,
							     ToyNumber);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}

std::unique_ptr<const RawBinnedDTYieldLikelihood>
BinnedDTData::GetDTYieldLikelihood(const std::string &Tag,
				   const Settings &settings,
				   int ToyNumber) const {
  if(!m_FullLikelihood) {
    return nullptr;
  }
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), Tag);
  const auto iter_Flavour = std::find(m_FlavourTags.begin(),
				      m_FlavourTags.end(), Tag);
  std::string TagCategory;
  if(iter_CP != m_CPTags.end()) {
    TagCategory = "CP";
  } else if(iter_SCMB != m_SCMBTags.end()) {
    TagCategory = "SCMB";
  } else if(iter_Flavour != m_FlavourTags.end()) {
    TagCategory = "Flavour";
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
  if(ToyNumber < 0) {
    std::cout << "Loading " << m_TagMode << " tag with full likelihood\n";
  }
  return std::make_unique<const RawBinnedDTYieldLikelihood>(Tag,
							    settings,
							    TagCategory,
							    ToyNumber);
}

std::unique_ptr<const BinnedDTYieldPrediction> BinnedDTData::GetDTPredictions(
  const std::string &Tag,
  const Settings &settings) const {
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), Tag);
  const auto iter_Flavour = std::find(m_FlavourTags.begin(),
				      m_FlavourTags.end(), Tag);
  if(iter_CP != m_CPTags.end()) {
    return std::make_unique<const BinnedCPTagYieldPrediction>(Tag, settings);
  } else if(iter_SCMB != m_SCMBTags.end()) {
    return std::make_unique<const BinnedSCMBTagYieldPrediction>(Tag, settings);
  } else if(iter_Flavour != m_FlavourTags.end()) {
    return std::make_unique<const BinnedFlavourTagYieldPrediction>(Tag, settings);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}

void BinnedDTData::SavePredictedBinYields(
  std::ofstream &File,
  const cisiFitterParameters &Parameters) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(Parameters);
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), m_TagMode);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), m_TagMode);
  const auto iter_Flavour = std::find(m_FlavourTags.begin(),
				      m_FlavourTags.end(), m_TagMode);
  if(iter_CP != m_CPTags.end()) {
    for(std::size_t Bin = 1; Bin <= Parameters.m_ci.size(); Bin++) {
      File << "DoubleTag_CP_KKpipi_vs_";
      File << m_TagMode << "_SignalBin" << Bin << " ";
      File << PredictedYields[Bin - 1] << "\n";
    }
  } else if(iter_SCMB != m_SCMBTags.end()) {
    // Assume K0pipi for now
    std::size_t Counter = 0;
    for(std::size_t TagBin = 1; TagBin <= 8; TagBin++) {
      for(int SignalBin = -Parameters.m_ci.size();
	  SignalBin <= static_cast<int>(Parameters.m_ci.size());
	  SignalBin++) {
	if(SignalBin == 0) {
	  continue;
	}
	File << "DoubleTag_SCMB_KKpipi_vs_";
	File << m_TagMode << "_SignalBin" << (SignalBin > 0 ? "P" : "M");
	File << TMath::Abs(SignalBin) << "_TagBin" << TagBin << " ";
	File << PredictedYields[Counter++] << "\n";
      }
    }
  } else if(iter_Flavour != m_FlavourTags.end()) {
    std::size_t Counter = 0;
    for(int Bin = -Parameters.m_ci.size();
	Bin <= static_cast<int>(Parameters.m_ci.size());
	Bin++) {
      if(Bin == 0) {
	continue;
      }
      File << "DoubleTag_Flavour_KKpipi_vs_";
      File << m_TagMode << "_SignalBin" << (Bin > 0 ? "P" : "M");
      File << TMath::Abs(Bin) << "_TagBin0 ";
      File << PredictedYields[Counter++] << "\n";
    }
  } else {
    throw std::runtime_error(m_TagMode + " is not a valid tag mode");
  }
  File << "\n";
}
