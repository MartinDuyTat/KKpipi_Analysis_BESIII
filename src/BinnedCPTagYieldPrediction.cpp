// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TMath.h"
#include"BinnedCPTagYieldPrediction.h"
#include"Settings.h"
#include"FlavourTags/KiCombiner.h"

BinnedCPTagYieldPrediction::BinnedCPTagYieldPrediction(const std::string &Tag,
						       const KiCombiner *Ki,
						       const Settings &settings):
  BinnedDTYieldPrediction(Tag, Ki, settings),
  m_FPlus(settings["FPlus_TagModes"].getD(Tag)) {
}

std::vector<double> BinnedCPTagYieldPrediction::GetPredictedBinYields(
  double BF_KKpipi,
  const std::vector<double> &ci,
  const std::vector<double> &si) const {
  const std::vector<double> KiVector = m_Ki->GetKi(ci, si);
  const std::size_t Size = KiVector.size()/2;
  TMatrixT<double> BinYields(Size, 1);
  for(std::size_t Bin = 1; Bin <= Size; Bin++) {
    const double Ki = KiVector[Bin - 1 + Size];
    const double Kbari = KiVector[-Bin + Size];
    const size_t i = Bin - 1;
    const double D0Yield = Ki;
    const double Dbar0Yield = Kbari;
    const double SqrtKK = TMath::Sqrt(Ki*Kbari);
    const double Interference = -2.0*SqrtKK*(2.0*m_FPlus - 1.0)*ci[i];
    const double UnnormalisedYield = D0Yield + Dbar0Yield + Interference;
    const double NormalisedYield = m_SingleTagYield*UnnormalisedYield*BF_KKpipi;
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
