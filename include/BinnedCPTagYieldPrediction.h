// Martin Duy Tat 4th October 2022
/**
 * BinnedCPTagYieldPrediction inherits from BinnedDTYieldPrediction
 * It calculates the yield prediction for CP tags
 */

#ifndef BINNEDCPTAGYIELDPREDICTION
#define BINNEDCPTAGYIELDPREDICTION

#include<vector>
#include<string>
#include"Settings.h"
#include"BinnedDTYieldPrediction.h"

class BinnedCPTagYieldPrediction: public BinnedDTYieldPrediction {
 public:
  /**
   * Constructor that stores the single tag yield, efficiency matrix and Ki
   * @param Tag The tag mode
   * @param settings The settings file
   */
  BinnedCPTagYieldPrediction(const std::string &Tag,
			     const Settings &settings);
  /**
   * Default virtual destructor
   */
  ~BinnedCPTagYieldPrediction() = default;
  /**
   * Function that returns the predicted bin yield
   */
  virtual std::vector<double> GetPredictedBinYields(
    double BF_KKpipi,
    const std::vector<double> &ci,
    const std::vector<double> &si,
    const std::vector<double> &Ki,
    const std::vector<double> &Kbari,
    double DeltaKpi) const override;
 private:
  /**
   * The CP-even fraction of the tag
   */
  const double m_FPlus;
};

#endif