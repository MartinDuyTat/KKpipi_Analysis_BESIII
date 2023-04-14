// Martin Duy Tat 13th February 2023
/**
 * BinnedSCMBTagYieldPrediction inherits from BinnedDTYieldPrediction
 * It calculates the yield prediction for SCMB tags
 */

#ifndef BINNEDSCMBTAGYIELDPREDICTION
#define BINNEDSCMBTAGYIELDPREDICTION

#include<vector>
#include<string>
#include"Settings.h"
#include"BinnedDTYieldPrediction.h"
#include"cisiK0pipi.h"

class BinnedSCMBTagYieldPrediction: public BinnedDTYieldPrediction {
 public:
  /**
   * Constructor that stores the single tag yield, efficiency matrix, signal mode Ki
   * It also loads Ki, ci and si for the tag mode
   * @param Tag The tag mode
   * @param settings The settings file
   */
  BinnedSCMBTagYieldPrediction(const std::string &Tag,
			       const Settings &settings);
  /**
   * Default virtual destructor
   */
  ~BinnedSCMBTagYieldPrediction() = default;
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
   * The tag mode
   */
  const std::string m_TagMode;
  /**
   * The Ki, ci and si of the tag mode
   */
  const cisiK0pipi m_Tag_Kicisi;
};

#endif
