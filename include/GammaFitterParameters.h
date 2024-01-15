// Martin Duy Tat 11th October 2023
/**
 * GammaFitterParameters is a simple struct that contains the fit parameters
 * for the \f$\gamma\f$ fit
 */

#ifndef GAMMAFITTERPARAMETERS
#define GAMMAFITTERPARAMETERS

#include"TMath.h"
#include"cisiFitterParameters.h"

struct GammaFitterParameters {
  /**
   * Constructor from a C-style array
   */
  GammaFitterParameters(const double *Parameters, std::size_t NBins):
    m_cisiParameters(Parameters, NBins),
    m_gamma(Parameters[4*NBins + 3]),
    //m_deltaB_dk(Parameters[4*NBins + 4]),
    //m_rB_dk(Parameters[4*NBins + 5]),
    m_rBcosDeltaB_dk(Parameters[4*NBins + 4]),
    m_rBsinDeltaB_dk(Parameters[4*NBins + 5]),
    m_deltaB_dk(TMath::ATan2(m_rBsinDeltaB_dk, m_rBcosDeltaB_dk)*180.0/TMath::Pi()),
    m_rB_dk(TMath::Sqrt(m_rBcosDeltaB_dk*m_rBcosDeltaB_dk +
			m_rBsinDeltaB_dk*m_rBsinDeltaB_dk)),
    //m_deltaB_dpi(Parameters[4*NBins + 6]),
    //m_rB_dpi(Parameters[4*NBins + 7]),
    m_xMinus(m_rB_dk*TMath::Cos((m_deltaB_dk - m_gamma)*TMath::Pi()/180.0)),
    m_yMinus(m_rB_dk*TMath::Sin((m_deltaB_dk - m_gamma)*TMath::Pi()/180.0)),
    m_xPlus(m_rB_dk*TMath::Cos((m_deltaB_dk + m_gamma)*TMath::Pi()/180.0)),
    m_yPlus(m_rB_dk*TMath::Sin((m_deltaB_dk + m_gamma)*TMath::Pi()/180.0)),
    //m_xXi((m_rB_dpi/m_rB_dk)*TMath::Cos((m_deltaB_dpi - m_deltaB_dk)*TMath::Pi()/180.0)),
    //m_yXi((m_rB_dpi/m_rB_dk)*TMath::Sin((m_deltaB_dpi - m_deltaB_dk)*TMath::Pi()/180.0)),
    m_xXi(Parameters[4*NBins + 6]),
    m_yXi(Parameters[4*NBins + 7]),
    m_Ri(std::vector<double>(Parameters + 4*NBins + 8, Parameters + 6*NBins + 7)),
    m_BMinusDKYield(Parameters[6*NBins + 7]),
    m_BPlusDKYield(Parameters[6*NBins + 8]),
    m_BMinusDpiYield(Parameters[6*NBins + 9]),
    m_BPlusDpiYield(Parameters[6*NBins + 10]) {
  }
  /**
   * The \f$c_i$ and \f$s_i\f$ parameters and other nuisance parameters
   */
  cisiFitterParameters m_cisiParameters;
  double m_gamma;
  double m_rBcosDeltaB_dk;
  double m_rBsinDeltaB_dk;
  double m_deltaB_dk;
  double m_rB_dk;
  double m_deltaB_dpi;
  double m_rB_dpi;
  /**
   * \f$x_-\f$
   */
  double m_xMinus;
  /**
   * \f$y_-\f$
   */
  double m_yMinus;
  /**
   * \f$x_+\f$
   */
  double m_xPlus;
  /**
   * \f$y_+\f$
   */
  double m_yPlus;
  /**
   * \f$x_\xi\f$
   */
  double m_xXi;
  /**
   * \f$y_\xi\f$
   */
  double m_yXi;
  /**
   * The \f$R_i\f$ parameters for LHCb
   */
  std::vector<double> m_Ri;
  /**
   * The \f$B^-\f$ yield for \f$DK\f$
   */
  double m_BMinusDKYield;
  /**
   * The \f$B^+\f$ yield for \f$DK\f$
   */
  double m_BPlusDKYield;
  /**
   * The \f$B^-\f$ yield for \f$D\pi\f$
   */
  double m_BMinusDpiYield;
  /**
   * The \f$B^+\f$ yield for \f$D\pi\f$
   */
  double m_BPlusDpiYield;
  /**
   * Function to update the value of \f$\gamma\f$
   */
  void UpdateGamma(double gamma) {
    m_gamma = gamma;
    m_xMinus = m_rB_dk*TMath::Cos((m_deltaB_dk - m_gamma)*TMath::Pi()/180.0);
    m_yMinus = m_rB_dk*TMath::Sin((m_deltaB_dk - m_gamma)*TMath::Pi()/180.0);
    m_xPlus = m_rB_dk*TMath::Cos((m_deltaB_dk + m_gamma)*TMath::Pi()/180.0);
    m_yPlus = m_rB_dk*TMath::Sin((m_deltaB_dk + m_gamma)*TMath::Pi()/180.0);
  }
};

#endif
