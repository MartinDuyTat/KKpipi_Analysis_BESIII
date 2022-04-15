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
#include"RooDataSet.h"
#include"RooMultiVarGaussian.h"
#include"Settings.h"
#include"cisiK0pipi.h"

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
     * Fit setup with all yields and efficiencies are stored in settings
     */
    Settings m_Settings;
    /**
     * Model predicted value of F+
     */
    const double m_FPlus_Model;
    /**
     * PDG value of KKpipi branching fraction
     */
    const double m_KKpipi_BF_PDG;
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
     * The KKpipi branching fraction for CP tags
     */
    RooRealVar m_KKpipi_BF_CP;
    /**
     * The KKpipi branching fraction for KSpipi tag
     */
    RooRealVar m_KKpipi_BF_KSpipi;
    /**
     * The KKpipi branching fraction for KLpipi tag
     */
    RooRealVar m_KKpipi_BF_KLpipi;
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
     * Perform a single fit to data
     */
    void DoSingleFit(RooMultiVarGaussian *Model);
    /**
     * Perform a single fit to a toy
     */
    void DoSingleToy(RooMultiVarGaussian *Model);
    /**
     * Perform many fits to data, potentially smearing stuff between each fit
     */
    //void DoManyFits(const RooMultiVarGaussian &Model, const RooDataSet &Data);
    /**
     * Perform many fits to toys
     */
    void DoManyToys(RooMultiVarGaussian *Model);
    /**
     * Vector of Gaussian constraint PDFs
     */
    RooArgSet m_GaussianConstraintPDFs;
    /**
     * Struct storing all the information about ci, si and Ki for KSpipi and KLpipi
     */
    cisiK0pipi m_cisi_K0pipi;
    /**
     * Resets the floating parameters to the F+ model value and PDG value of the BF
     */
    void ResetParameters();
};

#endif
