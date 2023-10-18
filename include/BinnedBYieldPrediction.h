// Martin Duy Tat 11th October 2023
/**
 * BinnedBYieldPrediction predicts the yields based on the value of \f$\gamma\f$
 * and other parameters
 */

#ifndef BINNEDBYIELDPREDICTION
#define BINNEDBYIELDPREDICTION

#include<vector>
#include<utility>
#include<array>
#include"Settings.h"
#include"GammaFitterParameters.h"

class BinnedBYieldPrediction {
 public:
  /**
   * Constructor that saves the number of bins
   */
  BinnedBYieldPrediction(const Settings &settings);
  /**
   * Function that predicts the \f$B^\pm\f$ yields
   * @param Parameters The fit parameters for \f$\gamma\f$
   */
  std::vector<double> GetPredictedYields(
    const GammaFitterParameters &Parameters) const;
 private:
  /**
   * Number of bins
   */
  const std::size_t m_NumberBins;
  /**
   * Flag that is true if the \f$c_i\f$ and \f$s_i$ are fixed (for debugging)
   */
  const bool m_cisiFixed;
  /**
   * Helper function to calculate bin yield for \f$B^+\f$ in positive bins
   */
  static double CalculateBinYield(double x, double y,
				  double Fi, double Fbari,
				  double ci, double si);
  /**
   * Helper function to pick out the correct \f$c_i\f$ and \f$s_i\f$
   */
  std::pair<double, double> Pickcisi(int Bin, int Charge,
				     const std::vector<double> &ci,
				     const std::vector<double> &si) const;
  /**
   * Helper function to pick out the correct \f$F_i\f$ and \f$\bar{F_i}\f$
   */
  static std::pair<double, double> PickFi(int Bin, int Charge,
					  const std::vector<double> &Fi,
					  const std::vector<double> &Fbari);
  /**
   * Helper function to convert \f$x_\xi\f$ and \f$y_\xi\f$ to \f$x_\pm\f$ and \f$y_\pm\f$
   */
  static std::pair<double, double> GetxyDpi(double x, double y,
					    double xXi, double yXi);
  /**
   * Helper function to get a set of bin yields
   */
  std::vector<double> GetSetOfBinYields(
    int Charge,
    double x,
    double y,
    const std::vector<double> &Fi,
    const std::vector<double> &Fbari,
    const std::vector<double> &ci,
    const std::vector<double> &si,
    double TotalYield) const;
  /**
   * Model predicted values of \f$c_i\f$
   */
  static constexpr std::array<double, 4> m_ci_Model{-0.362085, 0.808165, 0.835997, -0.380409};
  /**
   * Model predicted values of \f$s_i\f$
   */
  static constexpr std::array<double, 4> m_si_Model{-0.658906, -0.326282, 0.336997, 0.629657};
};

#endif
