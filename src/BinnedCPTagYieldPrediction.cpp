// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TMath.h"
#include"BinnedCPTagYieldPrediction.h"
#include"Settings.h"
#include"Utilities.h"
#include"cisiFitterParameters.h"

BinnedCPTagYieldPrediction::BinnedCPTagYieldPrediction(const std::string &Tag,
						       const Settings &settings):
  BinnedDTYieldPrediction(Tag, settings),
  m_FPlus(settings["FPlus_TagModes"].getD(Tag)) {
}

std::vector<double> BinnedCPTagYieldPrediction::GetPredictedBinYields(
  const cisiFitterParameters &Parameters) const {
  std::vector<double> Ki, Kbari;
  Utilities::ConvertRiToKi(Parameters.m_Ri, Ki, Kbari);
  const std::size_t Size = Parameters.m_ci.size();
  TMatrixT<double> BinYields(Size, 1);
  for(std::size_t Bin = 1; Bin <= Size; Bin++) {
    const size_t i = Bin - 1;
    const double D0Yield = Ki[i];
    const double Dbar0Yield = Kbari[i];
    const double SqrtKK = TMath::Sqrt(Ki[i]*Kbari[i]);
    const double Interference = -2.0*SqrtKK*(2.0*m_FPlus - 1.0)*Parameters.m_ci[i];
    const double UnnormalisedYield = D0Yield + Dbar0Yield + Interference;
    const double STYield = m_SingleTagYield/(1.0 - (2.0*m_FPlus - 1.0)*m_y);
    const double NormalisedYield =
      STYield*UnnormalisedYield*Parameters.m_BF_KKpipi;
    BinYields(i, 0) = NormalisedYield;
  }
  const auto EffMatrix = (1.0 - m_FPlus)*m_EfficiencyMatrix_CPEven
                       + m_FPlus*m_EfficiencyMatrix_CPOdd;
  const auto EffCorrBinYields = EffMatrix*BinYields;
  std::vector<double> FinalBinYields(Size);
  for(std::size_t i = 0; i < Size; i++) {
    FinalBinYields[i] = EffCorrBinYields(i, 0);
  }
  return FinalBinYields;
}
