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
   * Do the minimisation, save results and plot contours
   */
  void Minimise() const;
  /**
   * Run toy studies and save results
   */
  void RunToy(int ToyNumber) const;
  /**
   * Save predicted yields to a file for each tag
   */
  void SavePredictedYields() const;
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
   * Flag that is true if we also fit \f$\delta_{K\pi}\f$
   */
  bool m_FitDeltaKpi;
  /**
   * Helper function to set up minimiser
   */
  void SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const;
  /**
   * Helper function to generator values for toy generation
   */
  cisiFitterParameters GetGeneratorValues() const;
  /**
   * Function for drawing the contours of ci vs si
   * @param Minimiser The minimiser
   * @param ci The fitted values of ci
   * @param si The fitted values of si
   */
  void Plot_cisi(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
		 const std::vector<double> &ci,
		 const std::vector<double> &si) const;
  /**
   * Function for drawing the \f$\delta_{K\pi}\f$ contour
   * @param Minimiser The minimiser
   */
  void Plot_DeltaKpi(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
		     double rDcosDeltaKpi,
		     double rDsinDeltaKpi) const;
  /**
   * Function that saves the fit results
   * @param Minimiser The minimiser
   * @param Filename The filename where results are saved
   */
  void SaveFitResults(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
		      const std::string &Filename) const;
};

#endif
