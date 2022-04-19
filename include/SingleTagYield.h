// Martin Duy Tat 4th April 2021
/** 
 * SingleTagYield is a class for fitting the single tag yield
 * The yield is fitted with a signal shape from an exclusive signal MC sample, convolved with a Gaussian, and an Argus PDF for the combinatorial background
 * Any peaking backgrounds are fixed by Gaussians
 */

#ifndef SINGLETAGYIELD
#define SINGLETAGYIELD

#include<string>
#include<vector>
#include<utility>
#include<map>
#include"TTree.h"
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooExtendPdf.h"
#include"RooFitResult.h"
#include"RooAddPdf.h"
#include"RooDataSet.h"
#include"RooArgSet.h"
#include"Settings.h"
#include"RooShapes/FitShape.h"

class SingleTagYield {
  public:
    /**
     * Constructor that takes in the TTree object with the data and a TTree with the MC signal shape
     * @param DataTree TTree with signal events
     * @param MCSignalTree TTree with MC signal shape for fitting the signal
     */
    SingleTagYield(TTree *DataTree, TTree *MCSignalTree, const Settings &settings);
    /**
     * Destructor that deletes all the PDFs and corresponding variables for peaking backgrounds that were allocated on the heap
     */
    ~SingleTagYield();
    /**
     * Function that performs the fit of the single tag yield in RooFit
     */
    void FitYield();
  private:
    /**
     * TTree with signal events
     */
    TTree *m_DataTree;
    /**
     * TTree with exclusive signal MC for fitting the signal shape
     */
    TTree *m_MCSignalTree;
    /**
     * Settings object with the fit configuration
     */
    Settings m_Settings;
    /**
     * Parameters in fit
     */
    std::map<std::string, RooRealVar*> m_Parameters;
    /**
     * List of all PDFs
     */
    RooArgList m_ModelPDFs;
    /**
     * List of all PDF yields
     */
    RooArgList m_ModelYields;
    /**
     * Vector containing all the peaking backgrounds
     */
    std::vector<FitShape*> m_PeakingBackgrounds;
    /**
     * The full fit model of \f$\Delta E\f$
     */
    RooAddPdf* m_FullModel = nullptr;
    /**
     * The \f$m_{\rm BC}\f$ variable
     */
    RooRealVar m_MBC;
    /**
     * Weighting to account for luminosity scale
     */
    RooRealVar m_LuminosityWeight;
    /**
     * The fit results
     */
    RooFitResult *m_Result = nullptr;
    /**
     * Initial parameters before fit
     */
    RooArgSet *m_InitialParameters;
    /**
     * Initialize the signal shape
     */
    void InitializeSignalShape();
    /**
     * Initialize the combinatorial background shape
     */
    void InitializeArgus();
    /**
     * Initialize any peaking backgrounds
     */
    void InitializePeakingBackgrounds();
    /**
     * Initialize full mass shape
     */
    void InitializeFitShape();
    /**
     * Plot the single tag MBC fit
     * @param Data The data to be plotted
     */
    void PlotSingleTagYield(const RooDataSet &Data) const;
    /**
     * Function that saves the fit parameters to a text file
     */
    void SaveFitParameters() const;
    /**
     * Calculate the single tag yield and the uncertainty inside the signal region
     */
    std::pair<double, double> CalculateSingleTagYield() const;
    /**
     * Run the sPlot background subtract based on the fitted model
     * The result is saved as a TTree in a ROOT file
     * @param Data The dataset we want to apply sPlot on
     */
    void sPlotReweight(RooDataSet &Data);
    /**
     * Smear the Argus end point by 0.5 MeV
     */
    void SmearArgusEndPoint();
    /**
     * Smear the peaking backgrounds by their uncertainties
     */
    void SmearPeakingBackgrounds();
};

#endif
