// Martin Duy Tat 14th December 2021
/**
 * FPlusFitter is a class where the user can add double tag yields, normalized by the single tag yields, and this class will put together a fit model to extract the CP even fraction \f$F_+\f$
 */

#ifndef FPLUSFITTER
#define FPLUSFITTER

#include<vector>
#include"RooDataSet.h"
#include"RooArgSet.h"
#include"RooArgList.h"
#include"RooFitResult.h"
#include"RooRealVar.h"
#include"RooAbsPdf.h"
#include"RooArgSet.h"
#include"Settings.h"

class FPlusFitter {
  public:
    /**
     * Constructor that saves the settings
     * @param settings Settings containing yields and efficiencies
     */
    FPlusFitter(const Settings &settings);
    /**
     * Add a tag mode to combination
     * @param TagMode Name of tag mode
     */
    void AddTag(const std::string &TagMode);
    /**
     * Initalize multidimensional Gaussian and the dataset, then fit to obtain CP even fraction
     */
    void InitializeAndFit();
    /**
     * Get the F+ factors and if there is an uncertainty, Gaussian constrain them (by default these are fixed)
     * @param TagMode Name of tag mode
     */
    RooRealVar* GetFPlusTag(const std::string &TagMode);
  private:
    /**
     * Set containing all the double tag yields, normalized by the correponding single tag yields
     */
    RooArgSet m_NormalizedYields;
    /**
     * List of predicted yields for each normalized double tag yield
     */
    RooArgList m_PredictedYields;
    /**
     * Uncertainties of each normalized double tag yield
     */
    std::vector<double> m_Uncertainties;
    /**
     * The CP even fraction, which is the only floating parameter in the fit
     */
    RooRealVar m_FPlus;
    /**
     * The KKpipi branching fraction
     */
    RooRealVar m_KKpipi_BF;
    /**
     * Fit setup with all yields and efficiencies are stored in settings
     */
    Settings m_Settings;
    /**
     * Add datapoint containing normalized yield of CP tag to dataset
     * @param TagMode Name of tag mode
     */
    void AddMeasurement_CP(const std::string &TagMode);
    /**
     * Add formula for prediction of normalized yield of CP tag
     * @param TagMode Name of tag mode
     */
    void AddPrediction_CP(const std::string &TagMode);
    /**
     * Add datapoint containing binned, normalized yield of Kshh tag to dataset
     */
    void AddMeasurement_KShh(const std::string &TagMode);
    /**
     * Add formula for prediction of normalized yield of KShh tag
     */
    void AddPrediction_KShh(const std::string &TagMode);
    /**
     * Save the fit results
     */
    void SaveFitResults(RooFitResult *Result) const;
    /**
     * Vector of Gaussian constraint PDFs
     */
    RooArgSet m_GaussianConstraintPDFs;
};

#endif
