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

#include"RooRealVar.h"
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
     * Fit setup with all yields and efficiencies are stored in settings
     */
    Settings m_Settings;
    /**
     * Add datapoint containing normalized yield to dataset
     * @param TagMode Name of tag mode
     */
    void AddMeasurement(const std::string &TagMode);
    /**
     * Add formula for prediction of normalized yield
     * @param TagMode Name of tag mode
     */
    void AddPrediction(const std::string &TagMode);
};

#endif