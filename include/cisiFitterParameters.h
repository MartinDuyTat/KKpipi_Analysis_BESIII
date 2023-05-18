/**
 * cisiFitterParameters is a simple struct that contains the fit parameters
 */

#ifndef CISIFITTERPARAMETERS
#define CISIFITTERPARAMETERS

#include<vector>

struct cisiFitterParameters {
  /**
   * Constructor from a C-style array
   */
  cisiFitterParameters(const double *Parameters, std::size_t NBins):
    m_BF_KKpipi(Parameters[0]),
    m_BF_KKpipi_KLpipi(Parameters[4*NBins]),
    m_ci(std::vector<double>(Parameters + 0*NBins + 1, Parameters + 1*NBins + 1)),
    m_si(std::vector<double>(Parameters + 1*NBins + 1, Parameters + 2*NBins + 1)),
    m_Ri(std::vector<double>(Parameters + 2*NBins + 1, Parameters + 4*NBins)),
    m_rDcosDeltaKpi(Parameters[4*NBins + 1]),
    m_rDsinDeltaKpi(Parameters[4*NBins + 2]) {
    m_Ri.push_back(1.0);
  }
  /**
   * The \f$D^0\to KK\pi\pi\f$ branching fraction
   */
  double m_BF_KKpipi;
  /**
   * The \f$D^0\to KK\pi\pi\f$ branching fraction for the \f$K_L\pi\pi\f$ tag
   */
  double m_BF_KKpipi_KLpipi;
  /**
   * The \f$c_i\f$ parameters
   */
  std::vector<double> m_ci;
  /**
   * The \f$s_i\f$ parameters
   */
  std::vector<double> m_si;
  /**
   * The \f$R_i\f$ parameters
   */
  std::vector<double> m_Ri;
  /**
   * \f$r_D\cos(\delta_{K\pi})\f$
   */
  double m_rDcosDeltaKpi;
  /**
   * \f$r_D\sin(\delta_{K\pi})\f$
   */
  double m_rDsinDeltaKpi;
};

#endif
