// Martin Duy Tat 9th October 2022

#include<vector>
#include<memory>
#include<algorithm>
#include<stdexcept>
#include<iostream>
#include<iomanip>
#include"Eigen/Dense"
#include"TMatrixT.h"
#include"TMath.h"
#include"TRandom.h"
#include"BinnedDTData.h"
#include"RawBinnedDTYields.h"
#include"RawBinnedCPTagYields.h"
#include"RawBinnedSCMBTagYields.h"
#include"BinnedDTYieldPrediction.h"
#include"BinnedCPTagYieldPrediction.h"
#include"BinnedSCMBTagYieldPrediction.h"
#include"FlavourTags/KiCombiner.h"

BinnedDTData::BinnedDTData(const std::string &Tag,
			   const KiCombiner *Ki,
			   const Settings &settings):
  m_TagMode(Tag),
  m_DTYields(GetRawDTYields(Tag, settings)),
  m_DTPredictions(GetDTPredictions(Tag, Ki, settings)),
  m_SymmetricUncertainties(settings.getB("SymmetricUncertainties")),
  m_EnvelopeConstant(1.0),
  m_DisplayToyEfficiency(settings.getB("DisplayToyEfficiency")) {
}

double BinnedDTData::GetLogLikelihood(double BF_KKpipi,
				      const std::vector<double> &ci,
				      const std::vector<double> &si) const {
  /*for(std::size_t i = 0; i < ci.size(); i++) {
    if(ci[i]*ci[i] + si[i]*si[i] > 1.0) {
      return 1.0e30;
    }
  }*/
  const auto MeasuredYields = m_DTYields->GetDoubleTagYields();
  return GetLogLikelihood(BF_KKpipi, ci, si, MeasuredYields);
}

double BinnedDTData::GetLogLikelihood(
  double BF_KKpipi,
  const std::vector<double> &ci,
  const std::vector<double> &si,
  const std::vector<AsymmetricUncertainty> &MeasuredYields) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(BF_KKpipi, ci, si);
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

void BinnedDTData::GenerateToyYields(double BF_KKpipi,
				     const std::vector<double> &ci,
				     const std::vector<double> &si,
				     std::size_t StatsMultiplier) const {
  m_ToyDTYields.clear();
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(BF_KKpipi, ci, si);
  std::size_t TotalGeneratedToys = 0;
  std::size_t Counter = 0;
  while(Counter < StatsMultiplier) {
    const auto ToyYields = m_DTYields->GetToyYields(PredictedYields,
						    m_SymmetricUncertainties);
    const double LogLikelihood = GetLogLikelihood(BF_KKpipi,
						  ci, si,
						  ToyYields.first);
    const double Probability = TMath::Exp(-0.5*LogLikelihood);
    const double GenProbability = ToyYields.second;
    const double ToyProbability = Probability/GenProbability;
    const double RandomNumber = gRandom->Uniform(0, 1);
    if(ToyProbability/m_EnvelopeConstant > RandomNumber) {
      m_ToyDTYields.push_back(ToyYields.first);
      Counter++;
    }
    m_EnvelopeConstant = std::max(ToyProbability, m_EnvelopeConstant);
    TotalGeneratedToys++;
  }
  if(m_DisplayToyEfficiency) {
    std::cout << m_TagMode << " efficiency: ";
    std::cout << static_cast<double>(StatsMultiplier)/TotalGeneratedToys;
    std::cout << "\n";
    std::cout << m_TagMode << " envelope constant: ";
    std::cout << m_EnvelopeConstant;
    std::cout << "\n";
  }
}

double BinnedDTData::GetToyLogLikelihood(double BF_KKpipi,
					 const std::vector<double> &ci,
					 const std::vector<double> &si) const {
  double LL = 0.0;
  for(const auto &ToyDTYield : m_ToyDTYields) {
    LL += GetLogLikelihood(BF_KKpipi, ci, si, ToyDTYield);
  }
  return LL;
}

void BinnedDTData::PrintComparison(double BF_KKpipi,
				   const std::vector<double> &ci,
				   const std::vector<double> &si) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(BF_KKpipi, ci, si);
  const auto MeasuredYields = m_DTYields->GetDoubleTagYields();
  std::cout << "Fitted and predicted yield comparison for " << m_TagMode << ":\n";
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
  const double Std2 = Yield.PlusUncertainty*Yield.NegativeUncertainty;
  const double StdDiff = Yield.NegativeUncertainty - Yield.PlusUncertainty;
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
  const Settings &settings) const {
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), Tag);
  if(iter_CP != m_CPTags.end()) {
    return std::make_unique<const RawBinnedCPTagYields>(Tag, settings);
  } else if(iter_SCMB != m_SCMBTags.end()) {
    return std::make_unique<const RawBinnedSCMBTagYields>(Tag, settings);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}

std::unique_ptr<const BinnedDTYieldPrediction> BinnedDTData::GetDTPredictions(
  const std::string &Tag,
  const KiCombiner *Ki,
  const Settings &settings) const {
  const auto iter_CP = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  const auto iter_SCMB = std::find(m_SCMBTags.begin(), m_SCMBTags.end(), Tag);
  if(iter_CP != m_CPTags.end()) {
    return std::make_unique<const BinnedCPTagYieldPrediction>(Tag,
							      Ki,
							      settings);
  } else if(iter_SCMB != m_SCMBTags.end()) {
    return std::make_unique<const BinnedSCMBTagYieldPrediction>(Tag,
								Ki,
								settings);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}
