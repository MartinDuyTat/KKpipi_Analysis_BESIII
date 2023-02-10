// Martin Duy Tat 13th Octboer 2022
/**
 * cisiLikelihood is a class that combines all the likelihood contributions from each tag and creates a final likelihood that can be minimised
 */

#ifndef CISILIKELIHOOD
#define CISILIKELIHOOD

#include<vector>
#include"Settings.h"
#include"BinnedDTData.h"

class cisiLikelihood {
 public:
  /**
   * Constructor that takes in the settings and sets up the likelihoods
   * @param settings The settings file
   */
  cisiLikelihood(const Settings &settings);
  /**
   * Delete the copy constructor
   */
  cisiLikelihood(const cisiLikelihood &Likelihood) = delete;
  /**
   * Function that returns the likelihood based on ci and si values
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The cosine of the strong phases
   * @param si The sine of the strong phases
   */
  double CalculateLogLikelihood(double BF_KKpipi,
				const std::vector<double> &ci,
				const std::vector<double> &si) const;
  /**
   * Generate toy datasets for all tags
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The cosine of the strong phases used to generate toy
   * @param si The sine of the strong phases used to generate toy
   */
  void GenerateToy(double BF_KKpipi,
		   const std::vector<double> &ci,
		   const std::vector<double> &si) const;
  /**
   * Function that returns the likelihood based on ci and si values with toy data
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The cosine of the strong phases
   * @param si The sine of the strong phases
   */
  double CalculateToyLogLikelihood(double BF_KKpipi,
				   const std::vector<double> &ci,
				   const std::vector<double> &si) const;
  /**
   * Function that prints a comparison between predicted and measured yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   */
  void PrintComparison(double BF_KKpipi,
		       const std::vector<double> &ci,
		       const std::vector<double> &si) const;
 private:
  /**
   * Vector of all the tags
   */
  const std::vector<BinnedDTData> m_TagData;
  /**
   * The Ki parameters
   */
  std::vector<double> m_Ki;
  /**
   * The Kbari parameters
   */
  std::vector<double> m_Kbari;
  /**
   * Helper function that sets up all the tags
   * @param settings The settings file
   */
  std::vector<BinnedDTData> SetupTags(const Settings &settings) const;
  /**
   * Helper function that loads Ki and Kbari
   */
  std::pair<std::vector<double>, std::vector<double>>
  GetKi(const Settings &settings) const;
};

#endif
