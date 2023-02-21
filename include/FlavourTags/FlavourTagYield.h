// Martin Duy Tat 15th February 2023
/**
 * FlavourTagYield is a class that stores a binned flavour tag yield
 */

#ifndef FLAVOURTAGYIELD
#define FLAVOURTAGYIELD

#include<string>
#include<utility>
#include<vector>
#include"uncertainties/ureal.hpp"
#include"TMatrixT.h"
#include"Settings.h"

class FlavourTagYield {
 public:
  /**
   * Constructor that loads all the double tag yields and the covariance
   */
  FlavourTagYield(const std::string &TagMode, const Settings &settings);
  /**
   * Get the value of Ki, after DCS corrections, normalised by the ST yield
   */
  uncertainties::udouble GetKi(int Bin,
			       const std::vector<double> &ci,
			       const std::vector<double> &si) const;
  /**
   * Print the DCS corrections
   */
  void PrintDCS(const std::vector<double> &ci,
		const std::vector<double> &si) const;
 private:
  /**
   * The flavour tag name
   */
  const std::string m_TagMode;
  /**
   * Flag that is true if DCS corrections are included
   */
  const bool m_IncludeDCSCorrections;
  /**
   * The raw double tag yields
   */
  const std::vector<uncertainties::udouble> m_DTFlavourYields;
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
   * The covariance matrix of the charm parameters
   */
  const TMatrixT<double> m_CharmCovMatrix;
  /**
   * The raw single tag yield
   */
  const std::pair<double, double> m_RawSingleTagYield;
  /**
   * The single tag efficiency
   */
  const std::pair<double, double> m_SingleTagEfficiency;
  /**
   * Helper function that gets the DCS and efficiency corrected flavour tag yield
   * @param Bin Bin number
   */
  uncertainties::udouble GetDTFlavourYield(int Bin) const;
  /**
   * Helper function that calculates the DCS correction using 3 iterations
   */
  double CalculateDCSCorrection(int Bin,
				const std::vector<double> &ci,
				const std::vector<double> &si) const;
  /**
   * Helper function that calculates the DCS correction using some specific Ki
   */
  double CalculateDCSCorrection(int Bin,
				double Ki,
				double Kbari,
				const std::vector<double> &ci,
				const std::vector<double> &si) const;
  /**
   * Helper function that loads the flavour tag yields
   */
  std::vector<uncertainties::udouble> LoadFlavourYields(
    const std::string &TagMode,
    const Settings &settings) const;
  /**
   * Helper function that loads the flavour tag yield covariance matrix
   */
  TMatrixT<double> LoadCovarianceMatrix(const std::string &TagMode,
					const Settings &settings) const;
  /**
   * Helper function to load efficiency matrix of the raw yields
   */
  TMatrixT<double> LoadEfficiencyMatrix(const std::string &TagMode,
					const Settings &settings) const;
  /**
   * Helper function that loads the charm parameters
   */
  std::pair<double, double> LoadCharmParameter(
    const std::string &TagMode,
    const Settings &settings,
    const std::string &ParameterName) const;
  /**
   * Helper function that loads the charm parameter covariance matrix
   */
  TMatrixT<double> LoadCharmCovMatrix(const std::string &TagMode,
				      const Settings &settings) const;
  /**
   * Helper function that loads the single tag yield
   */
  std::pair<double, double> GetSTYield(const std::string &Tag,
				       const Settings &settings) const;
  /**
   * Helper function that loads the single tag efficiency
   */
  std::pair<double, double> GetSTEfficiency(const std::string &Tag,
					    const Settings &settings) const;
};

#endif
