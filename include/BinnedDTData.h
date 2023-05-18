// Martin Duy Tat 9th October 2022
/**
 * BinnedDTData stores derived RawBinnedDTYields and BinnedDTYieldPrediction objects for a specific tag, and from this it can compute the likelihood of this tag
 */

#ifndef BINNEDDOUBLETAGDATA
#define BINNEDDOUBLETAGDATA

#include<vector>
#include<memory>
#include<string>
#include<array>
#include<string_view>
#include<fstream>
#include"Eigen/Dense"
#include"TMatrixT.h"
#include"RawBinnedDTYields.h"
#include"BinnedDTYieldPrediction.h"
#include"RawBinnedDTYieldLikelihood.h"
#include"cisiFitterParameters.h"

class BinnedDTData {
 public:
  /**
   * Constructor that stores pointers to the raw yields and predicted yields
   * @param Tag The tag mode
   * @param Kbari The Kbari parameters normalised by the ST yield
   * @param settings The settings file
   */
  BinnedDTData(const std::string &Tag,
	       const Settings &settings);
  /**
   * Delete copy constructor
   */
  BinnedDTData(const BinnedDTData &binnedDTData) = delete;
  /**
   * Default move constructor
   */
  BinnedDTData(BinnedDTData &&binnedDTData) = default;
  /**
   * Calculate the log likelihood from this tag
   * @param Parameters The fit parameters
   */
  double GetLogLikelihood(const cisiFitterParameters &Parameters) const;
  /**
   * Calculate the log likelihood from this tag, using alternative yields
   * @param Parameter The fit parameters
   * @param MeasuredYields The measured yields
   */
  double GetLogLikelihood(
    const cisiFitterParameters &Parameters,
    const std::vector<AsymmetricUncertainty> &MeasuredYields) const;
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
   * Dump the predicted yields of this tag in the file
   */
  void SavePredictedBinYields(std::ofstream &File,
			      const cisiFitterParameters &Parameters) const;
 private:
  /**
   * The name of this tag mode
   */
  const std::string m_TagMode;
  /**
   * Flag that is true if the full likelihood of the yields is used
   */
  const bool m_FullLikelihood;
  /**
   * The measured raw double tag yields
   */
  mutable std::unique_ptr<const RawBinnedDTYields> m_DTYields;
  /**
   * The full likelihood of the double tag yields
   */
  mutable std::unique_ptr<const RawBinnedDTYieldLikelihood> m_DTYieldLikelihood;
  /**
   * The object that predicts double tag yields
   */
  std::unique_ptr<const BinnedDTYieldPrediction> m_DTPredictions;
  /**
   * The generated toy double tag yields
   */
  mutable std::vector<std::vector<AsymmetricUncertainty>> m_ToyDTYields;
  /**
   * Flag that is true if the uncertainties of the yields are symmetrized
   */
  const bool m_SymmetricUncertainties;
  /**
   * Reference to the settings
   */
  const Settings &m_Settings;
  /**
   * Helper function that calculates the standard deviation from asymmetric uncertainties
   */
  double GetAsymmetricStd(const AsymmetricUncertainty &Yield,
			  double Prediction) const;
  /**
   * Helper function that creates the inverse of the covariance matrix
   */
  Eigen::MatrixXd CreateInvCovarianceMatrix(
    const std::vector<double> &PredictedYields,
    const std::vector<AsymmetricUncertainty> &MeasuredYields,
    const TMatrixT<double> &CorrelationMatrix) const;
  /**
   * Helper function that creates the correct raw yield object
   * @param Tag The tag mode
   * @param settings The settings file
   * @param ToyNumber For toys, this labels toy number, otherwise set to -1
   */
  std::unique_ptr<const RawBinnedDTYields> GetRawDTYields(
    const std::string &Tag,
    const Settings &settings,
    int ToyNumber = -1) const;
  /**
   * Helper function that loads the full likelihood object
   * @param Tag The tag mode
   * @param settings The settings file
   * @param ToyNumber For toys, this labels toy number, otherwise set to -1
   */
  std::unique_ptr<const RawBinnedDTYieldLikelihood> GetDTYieldLikelihood(
    const std::string &Tag,
    const Settings &settings,
    int ToyNumber = -1) const;
  /**
   * Helper function that creates the correct yield prediction object
   * @param Tag The tag mode
   * @param STYieldFlavour The flavour single tag yield
   * @param Kbari The Kbari parameters
   * @param settings The settings file
   */
  std::unique_ptr<const BinnedDTYieldPrediction> GetDTPredictions(
    const std::string &Tag,
    const Settings &settings) const;
  /**
   * List of the CP tags used in the analysis
   */
  static constexpr std::array<std::string_view, 13> m_CPTags{{
    "KK", "KKPartReco", "pipi", "pipipi0", "pipipi0PartReco", "KSpi0pi0", "KLpi0",
    "KSpi0", "KSpi0PartReco", "KSeta", "KSetaPrimepipieta", "KSetaPrimerhogamma", "KSpipipi0"}};
  /**
   * List of the SCMB tags used in the analysis
   */
  static constexpr std::array<std::string_view, 3> m_SCMBTags{{
    "KSpipi", "KSpipiPartReco", "KLpipi"}};
  /**
   * List of the flavour tags used in the analysis
   */
  static constexpr std::array<std::string_view, 4> m_FlavourTags{{
    "Kpi", "Kpipi0", "Kpipipi", "KeNu"}};
};

#endif
