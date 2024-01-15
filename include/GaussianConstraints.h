// Martin Duy Tat 15th January 2024
/**
 * GaussianConstraints is a namespace containing all the Gaussian constraints
 */

#ifndef GAUSSIANCONSTRAINTS
#define GAUSSIANCONSTRAINTS

#include"TMatrixT.h"
#include"cisiFitterParameters.h"

namespace GaussianConstraints {
  /**
   * Gaussian constraint for \f$\delta_D^{K\pi}\f$
   */
  double AddDeltaKpiConstraint(const cisiFitterParameters &Parameters) {
    TMatrixT<double> InvCov(2, 2);
    InvCov(0, 0) = 10920.207;
    InvCov(1, 1) = 4952.4760;
    InvCov(0, 1) = -147.08102;
    InvCov(1, 0) = -147.08102;
    const double CosDeltaKpi_diff = Parameters.m_rDcosDeltaKpi + 0.0562;
    const double SinDeltaKpi_diff = Parameters.m_rDsinDeltaKpi + 0.011;
    const double CosDeltaKpi_diff2 = CosDeltaKpi_diff*CosDeltaKpi_diff;
    const double SinDeltaKpi_diff2 = SinDeltaKpi_diff*SinDeltaKpi_diff;
    const double CosSinDeltaKpi_diff = CosDeltaKpi_diff*SinDeltaKpi_diff;
    return CosDeltaKpi_diff2*InvCov(0, 0) +
           SinDeltaKpi_diff2*InvCov(1, 1) +
           CosSinDeltaKpi_diff*2*InvCov(0, 1);
  }
  /**
   * Calculate all the constraints
   * @param Parameters the fit parameters
   */
  double AddGaussianConstraints(const cisiFitterParameters &Parameters) {
    double LL = 0.0;
    LL += AddDeltaKpiConstraint(Parameters);
    return LL;
  }
}

#endif
