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
#include"TMatrixT.h"
#include"RawBinnedDTYields.h"
#include"BinnedDTYieldPrediction.h"

class BinnedDTData {
 public:
  /**
   * Constructor that stores pointers to the raw yields and predicted yields
   * @param Tag The tag mode
   * @param Ki The Ki parameters normalised by the ST yield
   * @param Kbari The Kbari parameters normalised by the ST yield
   * @param settings The settings file
   */
  BinnedDTData(const std::string &Tag,
	       const std::vector<double> &Ki,
	       const std::vector<double> &Kbari,
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
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   */
  double GetLogLikelihood(double BF_KKpipi,
			  const std::vector<double> &ci,
			  const std::vector<double> &si) const;
  /**
   * Calculate the log likelihood from this tag, using alternative yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   * @param MeasuredYields The measured yields
   */
  double GetLogLikelihood(
    double BF_KKpipi,
    const std::vector<double> &ci,
    const std::vector<double> &si,
    const std::vector<AsymmetricUncertainty> &MeasuredYields) const;
  /**
   * Generate toy yields using hit and miss method
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters used in toy generation
   * @param si The si parameters used in toy generation
   * @param StatsMultiplier The statistics multiplier
   */
  void GenerateToyYields(double BF_KKpipi,
			 const std::vector<double> &ci,
			 const std::vector<double> &si,
			 std::size_t StatsMultiplier) const;
  /**
   * Calculate the log likelihood from this tag using the generated toy yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   */
  double GetToyLogLikelihood(double BF_KKpipi,
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
   * The name of this tag mode
   */
  const std::string m_TagMode;
  /**
   * The measured raw double tag yields
   */
  std::unique_ptr<const RawBinnedDTYields> m_DTYields;
  /**
   * The object that predicts double tag yields
   */
  std::unique_ptr<const BinnedDTYieldPrediction> m_DTPredictions;
  /**
   * The generated toy double tag yields
   */
  mutable std::vector<std::vector<AsymmetricUncertainty>> m_ToyDTYields;
  /**
   * Helper function that calculates the standard deviation from asymmetric uncertainties
   */
  double GetAsymmetricStd(const AsymmetricUncertainty &Yield,
			  double Prediction) const;
  /**
   * Helper function that creates the covariance matrix
   */
  TMatrixT<double> CreateCovarianceMatrix(
    const std::vector<double> &PredictedYields,
    const std::vector<AsymmetricUncertainty> &MeasuredYields,
    TMatrixT<double> CorrelationMatrix) const;
  /**
   * Helper function that creates the correct raw yield object
   * @param Tag The tag mode
   * @param settings The settings file
   */
  std::unique_ptr<const RawBinnedDTYields> GetRawDTYields(
    const std::string &Tag,
    const Settings &settings) const;
  /**
   * Helper function that creates the correct yield prediction object
   * @param Tag The tag mode
   * @param STYieldFlavour The flavour single tag yield
   * @param Ki The Ki parameters
   * @param Kbari The Kbari parameters
   * @param settings The settings file
   */
  std::unique_ptr<const BinnedDTYieldPrediction> GetDTPredictions(
    const std::string &Tag,
    const std::vector<double> &Ki,
    const std::vector<double> &Kbari,
    const Settings &settings) const;
  /**
   * List of the CP tags used in the analysis
   */
  static constexpr std::array<std::string_view, 12> m_CPTags{{
    "KK", "KKPartReco", "pipi", "pipipi0", "KSpi0pi0", "KLpi0",
    "KSpi0", "KSpi0PartReco", "KSeta", "KSetaPrimepipieta", "KSetaPrimerhogamma", "KSomega"}};
  /**
   * List of the SCMB tags used in the analysis
   */
  static constexpr std::array<std::string_view, 2> m_SCMBTags{{
    "KSpipi", "KLpipi"}};
};

#endif
