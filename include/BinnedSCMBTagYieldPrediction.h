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
#include"cisiFitterParameters.h"

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
    const cisiFitterParameters &Parameters) const override;
 private:
  /**
   * The tag mode
   */
  const std::string m_TagMode;
  /**
   * The Ki, ci and si of the tag mode
   */
  cisiK0pipi* m_Tag_Kicisi;
  /**
   * The CP-even fraction of the tag
   */
  const double m_FPlus;
};

#endif
