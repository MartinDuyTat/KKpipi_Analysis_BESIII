// Martin Duy Tat 11th October 2023
/**
 * GammaLikelihood is a class that combines the \f$c_i\f$ and \f$s_i\f$ fit with
 * the LHCb yields to obtain a likelihood for \f$\gamma\f$
 */

#ifndef GAMMALIKELIHOOD
#define GAMMALIKELIHOOD

#include<vector>
#include<utility>
#include"TMatrixT.h"
#include"Settings.h"
#include"GammaFitterParameters.h"
#include"cisiLikelihood.h"
#include"BinnedBYieldPrediction.h"

class GammaLikelihood {
 public:
  /**
   * Constructor that sets up the settings and likelihoods
   * @param settings The settings file
   * @param Likelihood The likelihood for ci and si
   */
  GammaLikelihood(const Settings &settings,
		  const cisiLikelihood &Likelihood);
  /**
   * Delete the copy constructor
   */
  GammaLikelihood(const GammaLikelihood &Likelihood) = delete;
  /**
   * Function that returns the likelihood for \f$\gamm\f$
   */
  double CalculateLogLikelihood(const GammaFitterParameters &Parameters) const;
 private:
  /**
   * Function for loading the (inverse of the) covariance matrix
   * @param Filename Filename of yields from LHCb
   * @param NBins Number of bins
   */
  static TMatrixT<double> GetCovMatrix(const std::string &Filename,
				       std::size_t NBins);
  /**
   * Function for loading the bin yields from LHCb
   * @param Filename Filename of yields from LHCb
   * @param NBins Number of bins
   */
  static std::vector<double> GetBinnedBYields(const std::string &Filename,
					      std::size_t NBins);
  /**
   * Number of bins
   */
  const std::size_t m_NumberBins;
  /**
   * The inverse of the covariance matrix
   */
  const TMatrixT<double> m_InvCovMatrix;
  /**
   * The bin yields
   */
  const std::vector<double> m_BinnedBYields;
  /**
   * The binned \f$B^\pm\f$ yield predictions
   */
  const BinnedBYieldPrediction m_BinnedBYieldPrediction;
  /**
   * The \f$c_i$ and \f$s_i\f$ likelihood
   */
  const cisiLikelihood& m_cisiLikelihood;
  /**
   * The settings object
   */
  const Settings& m_Settings;
  /**
   * Flag that is true if the \f$c_i\f$ and \f$s_i$ are fixed (for debugging)
   */
  const bool m_cisiFixed;
  /**
   * Flag that is true if the \f$c_i\f$ and \f$s_i\f$ are approximated
   * with asymmetric uncertainties
   */
  const bool m_cisiAsymmetric;
  /**
   * The fitted values of \f$c_i\f$ for the \f$\gamma\f$ fit
   */
  const std::vector<double> m_ci_fitted;
  /**
   * The fitted uncertainties of \f$c_i\f$ for the \f$\gamma\f$ fit
   */
  const std::vector<double> m_ci_fitted_err;
  /**
   * The fitted values of \f$s_i\f$ for the \f$\gamma\f$ fit
   */
  const std::vector<double> m_si_fitted;
  /**
   * The fitted uncertainties of \f$s_i\f$ for the \f$\gamma\f$ fit
   */
  const std::vector<std::pair<double, double>> m_si_fitted_err;
  /**
   * Helper function that calculates the likelihood of \f$c_i\f$ and \f$s_i$
   */
  double GetcisiLikelihood(const cisiFitterParameters &Parameters) const;
  /**
   * Helper function that calculates the likelihood of \f$c_i\f$
   */
  double GetciLikelihood(double ci, std::size_t Parameter) const;
  /**
   * Helper function that calculates the likelihood of \f$s_i\f$
   */
  double GetsiLikelihood(double si, std::size_t Parameter) const;
};

#endif
