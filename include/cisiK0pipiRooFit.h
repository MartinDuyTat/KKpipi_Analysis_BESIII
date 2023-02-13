// Martin Duy Tat 14th April 2022
/**
 * cisiK0pipiRooFit is a struct for storing the RooFit variables of ci and si for K0pipi
 * By default all variables are set to be constant, but they are initalised with a multidimensional Gaussian to account for their uncertainties and correlations, and if set non-constant they will be floated but Gaussian constrained
 */

#include<vector>
#include<string>
#include"TMatrixT.h"
#include"RooAbsPdf.h"
#include"RooArgList.h"
#include"Settings.h"

struct cisiK0pipiRooFit {
  /**
   * Flag to check if parameters are initialised
   */
  bool Initialised = false;
  /**
   * Initialise ci, si and Ki
   * @param settings Settings containing paths to parameters
   */
  void Initialise(const Settings &settings);
  /**
   * Initialise Ki
   * @param settings Settings containing paths to parameters
   * @param TagMode Tag mode
   */
  void InitialiseKi(const Settings &settings, const std::string &TagMode);
  /**
   * Initialise ci and si
   * @param settings Settings containing paths to parameters
   */
  void Initialisecisi(const Settings &settings);
  /**
   * Parse correlation matrix
   * @param settings Settings containing paths to parameters
   */
  TMatrixT<double> ParseCorrelationMatrix(const Settings &settings) const;
  /**
   * ci and si parameters
   */
  RooArgList m_cisi;
  /**
   * Ki parameters for KSpipi
   */
  RooArgList m_Ki_KSpipi;
  /**
   * Ki parameters for KLpipi
   */
  RooArgList m_Ki_KLpipi;
  /**
   * PDFs of Gaussian constraints
   */
  std::vector<RooAbsPdf*> m_GaussianConstraintPDFs;
};
