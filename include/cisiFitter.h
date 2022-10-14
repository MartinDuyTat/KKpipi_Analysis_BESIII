// Martin Duy Tat 13th October 2022
/**
 * cisiFitter fits ci and si by minimising the likelihood
 */

#ifndef CISIFITTER
#define CISIFITTER

#include<vector>
#include"Settings.h"
#include"cisiLikelihood.h"

class cisiFitter {
 public:
  /**
   * Constructor that takes in the settings and sets up the likelihoods
   * @param settings The settings file
   */
  cisiFitter(const Settings &settings);
  /**
   * Delete the copy constructor
   */
  cisiFitter(const cisiFitter &fitter) = delete;
  /**
   * Do the minimisation
   */
  void Minimise() const;
 private:
  /**
   * The ci and si likelihood
   */
  const cisiLikelihood m_cisiLikelihood;
  /**
   * Number of bins in the binning scheme
   */
  const int m_NumberBins;
};

#endif
