// Martin Duy Tat 9th October 2022

#include<vector>
#include<memory>
#include<algorithm>
#include<stdexcept>
#include<iostream>
#include<iomanip>
#include"TMatrixT.h"
#include"TMath.h"
#include"TRandom.h"
#include"BinnedDTData.h"
#include"RawBinnedDTYields.h"
#include"RawBinnedCPTagYields.h"
#include"BinnedDTYieldPrediction.h"
#include"BinnedCPTagYieldPrediction.h"

BinnedDTData::BinnedDTData(const std::string &Tag,
			   const std::vector<double> &Ki,
			   const std::vector<double> &Kbari,
			   const Settings &settings):
  m_DTYields(GetRawDTYields(Tag, settings)),
  m_DTPredictions(GetDTPredictions(Tag, Ki, Kbari, settings)) {
}

double BinnedDTData::GetLogLikelihood(double BF_KKpipi,
				      const std::vector<double> &ci,
				      const std::vector<double> &si) const {
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
  auto CovMatrix = CreateCovarianceMatrix(PredictedYields,
					  MeasuredYields,
					  CorrelationMatrix);
  CovMatrix.Invert();
  double LogLikelihood = 0.0;
  std::size_t Size = ci.size();
  for(std::size_t i = 0; i < Size; i++) {
    const double Diff_i = MeasuredYields[i].Value - PredictedYields[i];
    for(std::size_t j = 0; j < Size; j++) {
      const double Diff_j = MeasuredYields[j].Value - PredictedYields[j];
      LogLikelihood += CovMatrix(i, j)*Diff_i*Diff_j;
    }
  }
  return LogLikelihood;
}

void BinnedDTData::GenerateToyYields(double BF_KKpipi,
				     const std::vector<double> &ci,
				     const std::vector<double> &si) const {
  while(true) {
    const auto ToyYields = m_DTYields->GetToyYields();
    const double LogLikelihood = GetLogLikelihood(BF_KKpipi, ci, si, ToyYields);
    const double Probability = TMath::Exp(-0.5*LogLikelihood);
    const double RandomNumber = gRandom->Uniform(0, 1);
    if(Probability > RandomNumber) {
      m_ToyDTYields = ToyYields;
      break;
    }
  }
}

double BinnedDTData::GetToyLogLikelihood(double BF_KKpipi,
					 const std::vector<double> &ci,
					 const std::vector<double> &si) const {
  return GetLogLikelihood(BF_KKpipi, ci, si, m_ToyDTYields);
}

void BinnedDTData::PrintComparison(double BF_KKpipi,
				   const std::vector<double> &ci,
				   const std::vector<double> &si) const {
  const auto PredictedYields =
    m_DTPredictions->GetPredictedBinYields(BF_KKpipi, ci, si);
  const auto MeasuredYields = m_DTYields->GetDoubleTagYields();
  std::cout << "Fitted and predicted yield comparison:\n";
  std::cout << std::left << std::setw(10) << "Fitted";
  std::cout << std::left << std::setw(10) << "Predicted\n";
  for(std::size_t i = 0; i < PredictedYields.size(); i++) {
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << MeasuredYields[i].Value;
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << PredictedYields[i] << "\n";
  }
}

double BinnedDTData::GetAsymmetricStd(const AsymmetricUncertainty &Yield,
				      double Prediction) const {
  const double Std2 = Yield.PlusUncertainty*Yield.NegativeUncertainty;
  const double StdDiff = Yield.PlusUncertainty - Yield.NegativeUncertainty;
  const double YieldDiff = Yield.Value - Prediction;
  const double Var = Std2 + StdDiff*YieldDiff;
  return TMath::Sqrt(Var);
}

TMatrixT<double> BinnedDTData::CreateCovarianceMatrix(
  const std::vector<double> &PredictedYields,
  const std::vector<AsymmetricUncertainty> &MeasuredYields,
  TMatrixT<double> CorrelationMatrix) const {
  std::size_t Size = PredictedYields.size();
  for(std::size_t i = 0; i < Size; i++) {
    const double Std_i = GetAsymmetricStd(MeasuredYields[i], PredictedYields[i]);
    for(std::size_t j = 0; j < Size; j++) {
      const double Std_j = GetAsymmetricStd(MeasuredYields[j], PredictedYields[j]);
      CorrelationMatrix(i, j) *= Std_i*Std_j;
    }
  }
  return CorrelationMatrix;
}

std::unique_ptr<const RawBinnedDTYields> BinnedDTData::GetRawDTYields(
    const std::string &Tag,
    const Settings &settings) const {
  const auto iter = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  if(iter != m_CPTags.end()) {
    return std::make_unique<const RawBinnedCPTagYields>(Tag, settings);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}

std::unique_ptr<const BinnedDTYieldPrediction> BinnedDTData::GetDTPredictions(
    const std::string &Tag,
    const std::vector<double> &Ki,
    const std::vector<double> &Kbari,
    const Settings &settings) const {
  const auto iter = std::find(m_CPTags.begin(), m_CPTags.end(), Tag);
  if(iter != m_CPTags.end()) {
    return std::make_unique<const BinnedCPTagYieldPrediction>(Tag,
							      Ki,
							      Kbari,
							      settings);
  } else {
    throw std::runtime_error(Tag + " is not a valid tag mode");
  }
}
