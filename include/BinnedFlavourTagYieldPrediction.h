// Martin Duy Tat 4th October 2022
/**
 * BinnedFlavourTagYieldPrediction inherits from BinnedDTYieldPrediction
 * It calculates the yield prediction for flavour tags
 */

#ifndef BINNEDFLAVOURTAGYIELDPREDICTION
#define BINNEDFLAVOURTAGYIELDPREDICTION

#include<vector>
#include<string>
#include"Settings.h"
#include"BinnedDTYieldPrediction.h"

class BinnedFlavourTagYieldPrediction: public BinnedDTYieldPrediction {
 public:
  /**
   * Constructor that stores the single tag yield, efficiency matrix and Ki
   * @param Tag The tag mode
   * @param settings The settings file
   */
  BinnedFlavourTagYieldPrediction(const std::string &Tag,
			     const Settings &settings);
  /**
   * Default virtual destructor
   */
  ~BinnedFlavourTagYieldPrediction() = default;
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
   * The ratio of DCS to CF amplitude, with uncertainty
   */
  const std::pair<double, double> m_rD;
  /**
   * The coherence factor, with uncertainty
   */
  const std::pair<double, double> m_R;
  /**
   * The charm strong phase difference between DCS and CF amplitudes
   */
  const std::pair<double, double> m_DeltaD;
  /**
   * Sine of the charm strong phase difference between DCS and CF amplitudes
   */
  const double m_SinDeltaD;
  /**
   * Cosine of the charm strong phase difference between DCS and CF amplitudes
   */
  const double m_CosDeltaD;
  /**
   * Flag that is true for the Kpi tag and if we float the strong phase
   */
  const bool m_FitDeltaKpi;
  /**
   * Helper function that loads the charm parameters
   */
  std::pair<double, double> LoadCharmParameter(
    const std::string &TagMode,
    const Settings &settings,
    const std::string &ParameterName) const;
};

#endif
