// Martin Duy Tat 13th October 2022
/**
 * cisiFitter fits ci and si by minimising the likelihood
 */

#ifndef CISIFITTER
#define CISIFITTER

#include<vector>
#include<string>
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
  const std::size_t m_NumberBins;
  /**
   * The settings
   */
  Settings m_Settings;
  /**
   * Helper function to set up minimiser
   */
  void SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const;
  /**
   * Helper function to load values of ci and si used for toy generation
   * @param c_or_s_or_K String indicating if we want ci or si or Ki or Kbari
   */
  std::vector<double> GetGeneratorcisi(const std::string &c_or_s_or_K) const;
  /**
   * Function for drawing the contours of ci vs si
   * @param Minimiser The minimiser
   * @param ci The fitted values of ci
   * @param si The fitted values of si
   */
  void Plot_cisi(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
		 const std::vector<double> &ci,
		 const std::vector<double> &si) const;
};

#endif
