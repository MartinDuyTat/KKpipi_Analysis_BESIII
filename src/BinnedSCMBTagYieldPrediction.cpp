// Martin Duy Tat 13th February 2023

#include<vector>
#include<string>
#include"TMath.h"
#include"BinnedSCMBTagYieldPrediction.h"
#include"Settings.h"
#include"FlavourTags/KiCombiner.h"

BinnedSCMBTagYieldPrediction::BinnedSCMBTagYieldPrediction(
  const std::string &Tag,
  const KiCombiner *Ki,
  const Settings &settings):
  BinnedDTYieldPrediction(Tag, Ki, settings),
  m_TagMode(Tag),
  m_Tag_Kicisi(settings) {
}

std::vector<double> BinnedSCMBTagYieldPrediction::GetPredictedBinYields(
  double BF_KKpipi,
  const std::vector<double> &ci,
  const std::vector<double> &si) const {
  const std::vector<double> KiVector = m_Ki->GetKi(ci, si);
  const std::size_t Size = KiVector.size()/2;
  const std::size_t TotalSize = 8*2*Size;
  TMatrixT<double> BinYields(TotalSize, 1);
  std::size_t BinYields_index = 0;
  for(std::size_t TagBin = 1; TagBin <= 8; TagBin++) {
    for(int SignalBin = -Size; SignalBin <= static_cast<int>(Size); SignalBin++) {
      if(SignalBin == 0) {
	continue;
      }
      const size_t j = TMath::Abs(SignalBin);
      const double Ki = KiVector[j - 1 + Size];
      const double Kbari = KiVector[-j + Size];
      const bool Conj = SignalBin < 0;
      // Get the Ki, ci and si information on the tag side
      const double TagKi = m_Tag_Kicisi.Get_Ki(TagBin, m_TagMode);
      const double TagKbari = m_Tag_Kicisi.Get_Kbari(TagBin, m_TagMode);
      const double Tagci = m_Tag_Kicisi.Get_ci(TagBin, m_TagMode);
      const double Tagsi = m_Tag_Kicisi.Get_si(TagBin, m_TagMode);
      // D0 yield
      const double D0Yield = TagKbari*(Conj ? Kbari : Ki);
      const double Dbar0Yield = TagKi*(Conj ? Ki : Kbari);
      const double SqrtKK = TMath::Sqrt(TagKi*TagKbari*Ki*Kbari);
      const double Interference = -2.0*SqrtKK*(Tagci*ci[j - 1] +
                                               Tagsi*si[j - 1]*(Conj ? -1 : +1));
      const double Sign = m_TagMode == "KLpipi" ? -1 : +1;
      const double UnnormalisedYield = D0Yield + Dbar0Yield + Sign*Interference;
      const double NormalisedYield = m_SingleTagYield*UnnormalisedYield*BF_KKpipi;
      BinYields(BinYields_index, 0) = NormalisedYield;
      BinYields_index++;
    }
  }
  /*const auto EffMatrix = (1.0 - m_FPlus)*m_EfficiencyMatrix_CPEven
                       + m_FPlus*m_EfficiencyMatrix_CPOdd;*/
  const auto EffMatrix = m_EfficiencyMatrix_CPEven;
  const auto EffCorrBinYields = EffMatrix*BinYields;
  std::vector<double> FinalBinYields(TotalSize);
  for(std::size_t i = 0; i < TotalSize; i++) {
    FinalBinYields[i] = EffCorrBinYields(i, 0);
  }
  return FinalBinYields;
}
