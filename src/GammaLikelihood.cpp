// Martin Duy Tat 11th October 2023

#include<vector>
#include<string_view>
#include"TMatrixT.h"
#include"TMatrixTSym.h"
#include"TMath.h"
#include"TFile.h"
#include"RooArgList.h"
#include"RooRealVar.h"
#include"GammaLikelihood.h"
#include"Settings.h"
#include"GammaFitterParameters.h"
#include"BinnedBYieldPrediction.h"

GammaLikelihood::GammaLikelihood(const Settings &settings,
				 const cisiLikelihood &Likelihood):
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")),
  m_InvCovMatrix(GetCovMatrix(settings.get("GammaResultsFile"), m_NumberBins)),
  m_BinnedBYields(GetBinnedBYields(settings.get("GammaResultsFile"), m_NumberBins)),
  m_BinnedBYieldPrediction(settings),
  m_cisiLikelihood(Likelihood),
  m_Settings(settings),
  m_cisiFixed(settings.getB("cisiFixed")),
  m_cisiAsymmetric(settings.getB("cisiAsymmetricUncertainties")),
  m_ci_fitted({settings.getD("c1_fitted"),
	       settings.getD("c2_fitted"),
	       settings.getD("c3_fitted"),
	       settings.getD("c4_fitted")}),
  m_ci_fitted_err({settings.getD("c1_fitted_err"),
	          settings.getD("c2_fitted_err"),
	          settings.getD("c3_fitted_err"),
	          settings.getD("c4_fitted_err")}),
  m_si_fitted({settings.getD("s1_fitted"),
	       settings.getD("s2_fitted"),
	       settings.getD("s3_fitted"),
	       settings.getD("s4_fitted")}),
  m_si_fitted_err({{settings.getD("s1_fitted_plus_err"),
	            settings.getD("s1_fitted_minus_err")},
	           {settings.getD("s2_fitted_plus_err"),
	            settings.getD("s2_fitted_minus_err")},
	           {settings.getD("s3_fitted_plus_err"),
	            settings.getD("s3_fitted_minus_err")},
	           {settings.getD("s4_fitted_plus_err"),
	            settings.getD("s4_fitted_minus_err")}}) {
}

double GammaLikelihood::CalculateLogLikelihood(
  const GammaFitterParameters &Parameters) const {
  const auto PredictedYields = m_BinnedBYieldPrediction.GetPredictedYields(Parameters);
  double LL = 0.0;
  for(std::size_t i = 0; i < 8*m_NumberBins; i++) {
    const double YieldDiff_i = m_BinnedBYields[i] - PredictedYields[i];
    for(std::size_t j = 0; j <= i; j++) {
      const double YieldDiff_j = m_BinnedBYields[j] - PredictedYields[j];
      const int SymmetryFactor = (i == j ? 1 : 2);
      LL += YieldDiff_i*YieldDiff_j*m_InvCovMatrix(i, j)*SymmetryFactor;
    }
  }
  if(!m_cisiFixed) {
    if(m_cisiAsymmetric) {
      LL += GetcisiLikelihood(Parameters.m_cisiParameters);
    } else {
      LL += m_cisiLikelihood.CalculateLogLikelihood(Parameters.m_cisiParameters);
    }
  }
  return LL;
}

TMatrixT<double> GammaLikelihood::GetCovMatrix(const std::string &Filename,
			      std::size_t NBins) {
  // Open file with fit results
  TFile File(Filename.c_str(), "READ");
  // Loop over all variables, save the order and uncertainties
  auto FloatingParameters = static_cast<RooArgList*>(File.Get("floating_param"));
  std::vector<std::size_t> Order;
  std::vector<double> Uncertainties;
  constexpr std::array<std::string_view, 2> BModes{"dk", "dpi"};
  constexpr std::array<std::string_view, 2> Charges{"minus", "plus"};
  for(const auto &BMode : BModes) {
    for(const auto &Charge : Charges) {
      for(int Bin = -static_cast<int>(NBins);
	  Bin <= static_cast<int>(NBins);
	  Bin++) {
	if(Bin == 0) {
	  continue;
	}
	std::string Name = "n_sig_" + std::string(BMode) + "_d2kkpipi_fiveL_";
	Name += std::string(Charge) + "_bin";
	Name += (Bin > 0 ? "p" : "m");
	Name += std::to_string(TMath::Abs(Bin));
	Order.push_back(FloatingParameters->index(Name.c_str()));
	auto YieldVar =
	  static_cast<RooRealVar*>(FloatingParameters->find(Name.c_str()));
	Uncertainties.push_back(YieldVar->getError());
      }
    }
  }
  // Double loop over all signal yields to save covariance matrix
  TMatrixTSym<double> *correlation_matrix = nullptr;
  File.GetObject("correlation_matrix", correlation_matrix);
  TMatrixT<double> CovMatrix(8*NBins, 8*NBins);
  std::size_t Index_i = 0;
  for(auto i : Order) {
    std::size_t Index_j = 0;
    for(auto j : Order) {
      CovMatrix(Index_i, Index_j) = (*correlation_matrix)(i, j);
      CovMatrix(Index_i, Index_j) *= Uncertainties[Index_i]*Uncertainties[Index_j];
      Index_j++;
    }
    Index_i++;
  }
  return CovMatrix.Invert();
}

std::vector<double> GammaLikelihood::GetBinnedBYields(const std::string &Filename,
				     std::size_t NBins) {
  std::vector<double> Yields;
  TFile File(Filename.c_str(), "READ");
  RooArgList *FloatingParameters = nullptr;
  File.GetObject("floating_param", FloatingParameters);
  constexpr std::array<std::string_view, 2> BModes{"dk", "dpi"};
  constexpr std::array<std::string_view, 2> Charges{"minus", "plus"};
  for(const auto &BMode : BModes) {
    for(const auto &Charge : Charges) {
      for(int Bin = -static_cast<int>(NBins);
	  Bin <= static_cast<int>(NBins);
	  Bin++) {
	if(Bin == 0) {
	  continue;
	}
	std::string Name = "n_sig_" + std::string(BMode) + "_d2kkpipi_fiveL_";
	Name += std::string(Charge) + "_bin";
	Name += (Bin > 0 ? "p" : "m");
	Name += std::to_string(TMath::Abs(Bin));
	auto YieldVar =
	  static_cast<RooRealVar*>(FloatingParameters->find(Name.c_str()));
	Yields.push_back(YieldVar->getVal());
      }
    }
  }
  return Yields;
}

double GammaLikelihood::GetcisiLikelihood(const cisiFitterParameters &Parameters) const {
  double LL = 0.0;
  for(std::size_t Parameter = 1; Parameter <= m_NumberBins; Parameter++) {
    LL += GetciLikelihood(Parameters.m_ci[Parameter - 1], Parameter);
    LL += GetsiLikelihood(Parameters.m_si[Parameter - 1], Parameter);
  }
  return LL;
}

double GammaLikelihood::GetciLikelihood(double ci, std::size_t Parameter) const {
  const double Diff = m_ci_fitted[Parameter - 1] - ci;
  return Diff*Diff/m_ci_fitted_err[Parameter - 1];
}

double GammaLikelihood::GetsiLikelihood(double si, std::size_t Parameter) const {
  const double Diff = m_si_fitted[Parameter - 1] - si;
  const double PlusUncertainty = m_si_fitted_err[Parameter - 1].first;
  const double MinusUncertainty = m_si_fitted_err[Parameter - 1].second;
  const double Std2 = PlusUncertainty*MinusUncertainty;
  const double StdDiff = MinusUncertainty - PlusUncertainty;
  const double Var = Std2 + StdDiff*Diff;
  return (Var > 0.0 ? Diff*Diff/Var : 1.0e10);
}
