// Martin Duy Tat 11th October 2023

#include<vector>
#include<utility>
#include<algorithm>
#include"TMath.h"
#include"BinnedBYieldPrediction.h"
#include"Settings.h"
#include"Utilities.h"
#include"GammaFitterParameters.h"

BinnedBYieldPrediction::BinnedBYieldPrediction(const Settings &settings):
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")),
  m_cisiFixed(settings.getB("cisiFixed")) {
}

std::vector<double> BinnedBYieldPrediction::GetPredictedYields(
  const GammaFitterParameters &Parameters) const {
  // Get the Fi parameters
  std::vector<double> Fi, Fbari;
  auto Ri = Parameters.m_Ri;
  Ri.push_back(1.0);
  Utilities::ConvertRiToKi(Ri, Fi, Fbari);
  // Start with B->DK with negative charge
  auto BMinusDKYields = GetSetOfBinYields(-1,
					  Parameters.m_xMinus,
					  Parameters.m_yMinus,
					  Fi,
					  Fbari,
					  Parameters.m_cisiParameters.m_ci,
					  Parameters.m_cisiParameters.m_si,
					  Parameters.m_BMinusDKYield);
  // Start with B->DK with positive charge
  auto BPlusDKYields = GetSetOfBinYields(+1,
					 Parameters.m_xPlus,
					 Parameters.m_yPlus,
					 Fi,
					 Fbari,
					 Parameters.m_cisiParameters.m_ci,
					 Parameters.m_cisiParameters.m_si,
					 Parameters.m_BPlusDKYield);
  // For B->Dpi we need to change xy parameterisation
  // Start with B->Dpi with negative charge
  auto xyMinus_Dpi = GetxyDpi(Parameters.m_xMinus,
			      Parameters.m_yMinus,
			      Parameters.m_xXi,
			      Parameters.m_yXi);
  auto BMinusDpiYields = GetSetOfBinYields(-1,
					   xyMinus_Dpi.first,
					   xyMinus_Dpi.second,
					   Fi,
					   Fbari,
					   Parameters.m_cisiParameters.m_ci,
					   Parameters.m_cisiParameters.m_si,
					   Parameters.m_BMinusDpiYield);
  // Start with B->Dpi with positive charge
  auto xyPlus_Dpi = GetxyDpi(Parameters.m_xPlus,
			     Parameters.m_yPlus,
			     Parameters.m_xXi,
			     Parameters.m_yXi);
  auto BPlusDpiYields = GetSetOfBinYields(+1,
					  xyPlus_Dpi.first,
					  xyPlus_Dpi.second,
					  Fi,
					  Fbari,
					  Parameters.m_cisiParameters.m_ci,
					  Parameters.m_cisiParameters.m_si,
					  Parameters.m_BPlusDpiYield);
  // Join the vectors and return
  BMinusDKYields.insert(BMinusDKYields.end(),
			BPlusDKYields.begin(),
			BPlusDKYields.end());
  BMinusDKYields.insert(BMinusDKYields.end(),
			BMinusDpiYields.begin(),
			BMinusDpiYields.end());
  BMinusDKYields.insert(BMinusDKYields.end(),
			BPlusDpiYields.begin(),
			BPlusDpiYields.end());
  return BMinusDKYields;
}

double BinnedBYieldPrediction::CalculateBinYield(double x, double y,
						 double Fi, double Fbari,
						 double ci, double si) {
  return Fbari + Fi*(x*x + y*y) + 2*TMath::Sqrt(Fi*Fbari)*(x*ci - y*si);
}

std::pair<double, double> BinnedBYieldPrediction::Pickcisi(
  int Bin, int Charge,
  const std::vector<double> &ci,
  const std::vector<double> &si) const {
  int AbsBin = TMath::Abs(Bin);
  double c = ci[AbsBin - 1];
  double s = si[AbsBin - 1];
  if(m_cisiFixed) {
    c = m_ci_Model[AbsBin - 1];
    s = m_si_Model[AbsBin - 1];
  }
  if((Charge > 0 && Bin < 0) || (Charge < 0 && Bin > 0)) {
    s = -s;
  }
  return std::make_pair(c, s);
}

std::pair<double, double> BinnedBYieldPrediction::PickFi(
  int Bin, int Charge,
  const std::vector<double> &Fi,
  const std::vector<double> &Fbari) {
  int AbsBin = TMath::Abs(Bin);
  double F = Fi[AbsBin - 1];
  double Fbar = Fbari[AbsBin - 1];
  if((Charge > 0 && Bin < 0) || (Charge < 0 && Bin > 0)) {
    std::swap(F, Fbar);
  }
  return std::make_pair(F, Fbar);
}

std::pair<double, double> BinnedBYieldPrediction::GetxyDpi(
  double x, double y, double xXi, double yXi) {
  double x_Dpi = xXi*x - yXi*y;
  double y_Dpi = xXi*y + yXi*x;
  return std::make_pair(x_Dpi, y_Dpi);
}
  
std::vector<double> BinnedBYieldPrediction::GetSetOfBinYields(
  int Charge,
  double x,
  double y,
  const std::vector<double> &Fi,
  const std::vector<double> &Fbari,
  const std::vector<double> &ci,
  const std::vector<double> &si,
  double TotalYield) const {
  std::vector<double> Yields;
  for(int Bin = -static_cast<int>(m_NumberBins);
      Bin <= static_cast<int>(m_NumberBins);
      Bin++) {
    if(Bin == 0) {
      continue;
    }
    auto cisi = Pickcisi(Bin, Charge, ci, si);
    auto FiFbari = PickFi(Bin, Charge, Fi, Fbari);
    Yields.push_back(CalculateBinYield(x, y,
				       FiFbari.first,
				       FiFbari.second,
				       cisi.first,
				       cisi.second));
  }
  double Sum = std::accumulate(Yields.begin(), Yields.end(), 0.0);
  std::transform(Yields.begin(),
		 Yields.end(),
		 Yields.begin(),
		 [=] (double a) { return a*TotalYield/Sum; });
  return Yields;
}
