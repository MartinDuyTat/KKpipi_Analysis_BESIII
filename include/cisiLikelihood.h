// Martin Duy Tat 13th Octboer 2022
/**
 * cisiLikelihood is a class that combines all the likelihood contributions from each tag and creates a final likelihood that can be minimised
 */

#ifndef CISILIKELIHOOD
#define CISILIKELIHOOD

#include<vector>
#include<fstream>
#include"Settings.h"
#include"BinnedDTData.h"
#include"cisiFitterParameters.h"

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
   * @param Parameters The fit parameters
   */
  double CalculateLogLikelihood(const cisiFitterParameters &Parameters) const;
  /**
   * Load toy dataset
   * @param ToyNumber The toy dataset number
   */
  void LoadToyDataset(int ToyNumber) const;
  /**
   * Function that prints a comparison between predicted and measured yields
   * @param Parameters The fit parameters
   */
  void PrintComparison(const cisiFitterParameters &Parameters) const;
  /**
   * Function that prints the Ki
   */
  void PrintFinalKi(const std::vector<double> &Ri) const;
  /**
   * Dump the predicted yields for each tag in the file
   */
  void SavePredictedBinYields(std::ofstream &File,
			      const cisiFitterParameters &Parameters) const;
 private:
  /**
   * Vector of all the tags
   */
  const std::vector<BinnedDTData> m_TagData;
  /**
   * Helper function that sets up all the tags
   * @param settings The settings file
   */
  std::vector<BinnedDTData> SetupTags(const Settings &settings) const;
};

#endif
