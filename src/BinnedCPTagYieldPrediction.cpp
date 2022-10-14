// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TMath.h"
#include"BinnedCPTagYieldPrediction.h"
#include"Settings.h"

BinnedCPTagYieldPrediction::BinnedCPTagYieldPrediction(const std::string &Tag,
						       const std::vector<double> &Ki,
						       const std::vector<double> &Kbari,
						       const Settings &settings):
  BinnedDTYieldPrediction(Tag, Ki, Kbari, settings),
  m_FPlus(settings["FPlus_TagModes"].getD(Tag)) {
}

std::vector<double> BinnedCPTagYieldPrediction::GetPredictedBinYields(
  const std::vector<double> &ci, const std::vector<double>&) const {
  const std::size_t Size = m_Ki.size();
  TMatrixT<double> BinYields(Size, 1);
  for(std::size_t Bin = 1; Bin <= Size; Bin++) {
    const size_t i = Bin - 1;
    const double D0Yield = m_Ki[i];
    const double Dbar0Yield = m_Kbari[i];
    const double SqrtKK = TMath::Sqrt(m_Ki[i]*m_Kbari[i]);
    const double Interference = -2.0*SqrtKK*(2.0*m_FPlus - 1.0)*ci[i];
    const double UnnormalisedYield = D0Yield + Dbar0Yield + Interference;
    const double NormalisedYield = m_SingleTagYield*UnnormalisedYield;
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
