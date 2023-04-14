// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TMath.h"
#include"BinnedFlavourTagYieldPrediction.h"
#include"Settings.h"
#include"Utilities.h"

BinnedFlavourTagYieldPrediction::BinnedFlavourTagYieldPrediction(const std::string &Tag,
						       const Settings &settings):
  BinnedDTYieldPrediction(Tag, settings),
  m_rD(LoadCharmParameter(Tag, settings, "rD")),
  m_R(LoadCharmParameter(Tag, settings, "R")),
  m_DeltaD(LoadCharmParameter(Tag, settings, "deltaD")),
  m_SinDeltaD(TMath::Sin(TMath::Pi()*m_DeltaD.first)/180.0),
  m_CosDeltaD(TMath::Cos(TMath::Pi()*m_DeltaD.first)/180.0),
  m_FitDeltaKpi(Tag == "Kpi" && settings.getB("FitDeltaKpi")) {
}

std::vector<double> BinnedFlavourTagYieldPrediction::GetPredictedBinYields(
  double BF_KKpipi,
  const std::vector<double> &ci,
  const std::vector<double> &si,
  const std::vector<double> &Ri,
  double DeltaKpi) const {
  std::vector<double> Ki, Kbari;
  Utilities::ConvertRiToKi(Ri, Ki, Kbari);
  const std::size_t Size = ci.size();
  TMatrixT<double> BinYields(2*Size, 1);
  double CosDeltaD, SinDeltaD;
  if(m_FitDeltaKpi) {
    CosDeltaD = TMath::Cos(TMath::Pi()*DeltaKpi/180.0);
    SinDeltaD = TMath::Sin(TMath::Pi()*DeltaKpi/180.0);
  } else {
    CosDeltaD = m_CosDeltaD;
    SinDeltaD = m_SinDeltaD;
  }
  for(int Bin = -Size; Bin <= static_cast<int>(Size); Bin++) {
    if(Bin == 0) {
      continue;
    }
    const size_t i = TMath::Abs(Bin) - 1;
    double K, Kbar;
    if(Bin > 0) {
      K = Ki[i];
      Kbar = Kbari[i];
    } else { 
      K = Kbari[i];
      Kbar = Ki[i];
    }
    const double D0Yield = Kbar;
    const double Dbar0Yield = K*m_rD.first*m_rD.first;
    const double SqrtKK = TMath::Sqrt(K*Kbar);
    const int Sign = Bin > 0 ? +1 : -1;
    const double PhaseTerm = ci[i]*CosDeltaD + Sign*si[i]*SinDeltaD;
    const double Interference = -2.0*m_rD.first*m_R.first*SqrtKK*PhaseTerm;
    const double UnnormalisedYield = D0Yield + Dbar0Yield + Interference;
    const double NormalisedYield = m_SingleTagYield*UnnormalisedYield*BF_KKpipi;
    if(Bin < 0) {
      BinYields(Bin + Size, 0) = NormalisedYield;
    } else {
      BinYields(Bin + Size - 1, 0) = NormalisedYield;
    }
  }
  const auto EffCorrBinYields = m_EfficiencyMatrix*BinYields;
  std::vector<double> FinalBinYields(2*Size);
  for(std::size_t i = 0; i < 2*Size; i++) {
    FinalBinYields[i] = EffCorrBinYields(i, 0);
  }
  return FinalBinYields;
}

std::pair<double, double>
BinnedFlavourTagYieldPrediction::LoadCharmParameter(
  const std::string &TagMode,
  const Settings &settings,
  const std::string &ParameterName) const {
  const std::string Name = TagMode + "_DCS_Parameters";
  const double Value = settings[Name].getD(ParameterName);
  const double Error = settings[Name].getD(ParameterName + "_err");
  return std::make_pair(Value, Error);
}
