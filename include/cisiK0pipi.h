// Martin Duy Tat 13th February 2023
/**
 * cisiK0pipi is a class for storing the Ki, ci and si for K0pipi
 * The variables can also be smeared, according to their covariance matrix
 */

#include<vector>
#include"TMatrixT.h"
#include"Settings.h"

class cisiK0pipi {
 public:
  /**
   * Get the object, if it doesn't exist create it
   */
  static cisiK0pipi* Initialise(const Settings &settings);
  /**
   * Get Ki in bin i
   */
  double Get_Ki(std::size_t Bin, const std::string &TagMode) const;
  /**
   * Get Kbari in bin i
   */
  double Get_Kbari(std::size_t Bin, const std::string &TagMode) const;
  /**
   * Get ci in bin i
   */
  double Get_ci(std::size_t Bin, const std::string &TagMode) const;
  /**
   * Get si in bin i
   */
  double Get_si(std::size_t Bin, const std::string &TagMode) const;
 private:
  /**
   * Constructor that initialises all the parameters
   * It's private to avoid multiple instances
   */
  cisiK0pipi(const Settings &settings);
  /**
   * Initialise Ki
   * @param settings Settings containing paths to parameters
   * @param TagMode Tag mode
   * @return Ki parameters and their uncertainties
   */
  std::vector<double> InitialiseKi(const Settings &settings,
						      const std::string &TagMode);
  /**
   * Initialise ci and si
   * @param settings Settings containing paths to parameters
   */
  std::vector<double> Initialisecisi(const Settings &settings);
  /**
   * Parse correlation matrix
   * @param settings Settings containing paths to parameters
   */
  TMatrixT<double> ParseCorrelationMatrix(const Settings &settings) const;
  /**
   * The object itself
   */
  static cisiK0pipi *m_Instance;
  /**
   * Ki parameters for KSpipi
   */
  const std::vector<double> m_Ki_KSpipi;
  /**
   * Ki parameters for KLpipi
   */
  const std::vector<double> m_Ki_KLpipi;
  /**
   * ci and si parameters
   */
  const std::vector<double> m_cisi;
};
