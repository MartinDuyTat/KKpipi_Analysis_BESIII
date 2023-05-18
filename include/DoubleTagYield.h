// Martin Duy Tat 26th November 2021
/**
 * DoubleTagYield is a class that takes in the TTree with the double tag events and performs a simultaneous fit in each bin to determine the signal yield
 */

#ifndef DOUBLETAGYIELD
#define DOUBLETAGYIELD

#include"TTree.h"
#include"RooRealVar.h"
#include"RooFitResult.h"
#include"BinnedDataLoader.h"
#include"BinnedFitModel.h"
#include"Settings.h"
#include"Category.h"

class DoubleTagYield {
  public:
    /**
     * Constructor that takes in the settings and double tag events
     * @param settings The fit settings
     * @param Tree TTree with the double tag events in data
     */
    DoubleTagYield(const Settings &settings, TTree *Tree);
    /**
     * Perform simultaneous fit to determine double tag yields
     */
    void DoFit();
  private:
    /**
     * The fit variable
     */
    RooRealVar m_SignalMBC;
    /**
     * The fit settings
     */
    Settings m_Settings;
    /**
     * TTree with double tag events
     */
    TTree *m_Tree;
    /**
     * The data to be fitted
     */
    BinnedDataLoader m_DataLoader;
    /**
     * The fit model
     */
    BinnedFitModel m_FitModel;
    /**
     * Save the full likelihood function
     */
    void SaveLikelihood(const std::string &Filaname, RooDataSet *DataSet);
    /**
     * Perform toy studies
     */
    void DoToyFits();
    /**
     * Perform systematics studies
     */
    void DoSystematicsFits();
    /**
     * Plot projections of each bin in the fit
     */
    void PlotProjections();
    /**
     * Save signal yields from fit
     */
    void SaveSignalYields(const std::string &Filename,
			  RooFitResult *Result) const;
    /**
     * Helper function to find sideband yield with correctly reconstructed signal side and incorrect tag side reconstruction
     * Only use for fully reconstructed tags
     */
    double GetSidebandYield(int SignalBin, int TagBin) const;
    /**
     * Function that performs the sPlot background subtraction
     * @param Data DataSet
     * @param FitModel Fit model
     */
    void sPlotReweight(RooDataSet &Data);
    /**
     * Initial parameters before fit
     */
    RooArgSet *m_InitialParameters;
    /**
     * Parameters after nominal fit
     */
    RooArgSet *m_ParametersAfterFit;
};

#endif
