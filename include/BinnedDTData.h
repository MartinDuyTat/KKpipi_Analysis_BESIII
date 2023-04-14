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
#include"Eigen/Dense"
#include"TMatrixT.h"
#include"RawBinnedDTYields.h"
#include"BinnedDTYieldPrediction.h"
#include"FlavourTags/KiCombiner.h"

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
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   * @param Ki The Ki parameters
   * @param Kbari The Kbari parameters
   * @param DeltaKpi The D0->Kpi strong phase
   */
  double GetLogLikelihood(double BF_KKpipi,
			  const std::vector<double> &ci,
			  const std::vector<double> &si,
			  const std::vector<double> &Ki,
			  const std::vector<double> &Kbari,
			  double DeltaKpi) const;
  /**
   * Calculate the log likelihood from this tag, using alternative yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   * @param Ki The Ki parameters
   * @param Kbari The Kbari parameters
   * @param DeltaKpi The D0->Kpi strong phase
   * @param MeasuredYields The measured yields
   */
  double GetLogLikelihood(
    double BF_KKpipi,
    const std::vector<double> &ci,
    const std::vector<double> &si,
    const std::vector<double> &Ki,
    const std::vector<double> &Kbari,
    double DeltaKpi,
    const std::vector<AsymmetricUncertainty> &MeasuredYields) const;
  /**
   * Generate toy yields using hit and miss method
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters used in toy generation
   * @param si The si parameters used in toy generation
   * @param Ki The Ki parameters used in toy generation
   * @param Kbari The Kbari parameters used in toy generation
   * @param DeltaKpi The D0->Kpi strong phase used in toy generation
   * @param StatsMultiplier The statistics multiplier
   */
  void GenerateToyYields(double BF_KKpipi,
			 const std::vector<double> &ci,
			 const std::vector<double> &si,
			 const std::vector<double> &Ki,
			 const std::vector<double> &Kbari,
			 double DeltaKpi,
			 std::size_t StatsMultiplier) const;
  /**
   * Calculate the log likelihood from this tag using the generated toy yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   * @param Ki The Ki parameters
   * @param Kbari The Kbari parameters
   * @param DeltaKpi The D0->Kpi strong phase
   */
  double GetToyLogLikelihood(double BF_KKpipi,
			     const std::vector<double> &ci,
			     const std::vector<double> &si,
			     const std::vector<double> &Ki,
			     const std::vector<double> &Kbari,
			     double DeltaKpi) const;
  /**
   * Function that prints a comparison between predicted and measured yields
   * @param BF_KKpipi The KKpipi branching fraction
   * @param ci The ci parameters
   * @param si The si parameters
   * @param Ki The Ki parameters
   * @param Kbari The Kbari parameters
   * @param DeltaKpi The D0->Kpi strong phase
   */
  void PrintComparison(double BF_KKpipi,
		       const std::vector<double> &ci,
		       const std::vector<double> &si,
		       const std::vector<double> &Ki,
		       const std::vector<double> &Kbari,
		       double DeltaKpi) const;
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
   * Flag that is true if the uncertainties of the yields are symmetrized
   */
  const bool m_SymmetricUncertainties;
  /**
   * The constant to ensure that the Gaussian envelope covers the whole PDF
   */
  mutable double m_EnvelopeConstant;
  /**
   * Flag that is true if the toy generation efficiency is displayed
   */
  const bool m_DisplayToyEfficiency;
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
   */
  std::unique_ptr<const RawBinnedDTYields> GetRawDTYields(
    const std::string &Tag,
    const Settings &settings) const;
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
    "KSpi0", "KSpi0PartReco", "KSeta", "KSetaPrimepipieta", "KSetaPrimerhogamma", "KSomega"}};
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
  /**
   * Helper function to find the general Poisson parameter for a non-integer yield
   * @param Yield The yield we want to find the Poisson parameter of
   */
  double FindPoissonParameter(double Yield) const;
  /**
   * Helper function to find the asymmetric Poisson uncertainties of a yield
   * @param Yield The yield we want to find the asymmetric uncertainties of
   */
  std::pair<double, double> GetAsymmetricUncertainties(double Yield) const;
};

#endif
