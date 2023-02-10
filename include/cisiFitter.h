// Martin Duy Tat 13th October 2022
/**
 * cisiFitter fits ci and si by minimising the likelihood
 */

#ifndef CISIFITTER
#define CISIFITTER

#include<vector>
#include"Minuit2/Minuit2Minimizer.h"
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
  /**
   * Run toy studies and save results to a TTree
   */
  void RunToys() const;
 private:
  /**
   * The ci and si likelihood
   */
  const cisiLikelihood m_cisiLikelihood;
  /**
   * Number of bins in the binning scheme
   */
  const int m_NumberBins;
  /**
   * The settings
   */
  Settings m_Settings;
  /**
   * Helper function to set up minimiser
   */
  void SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const;
};

#endif
