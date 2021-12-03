// Martin Duy Tat 26th November 2021
/**
 * DoubleTagYield is a class that takes in the TTree with the double tag events and performs a simultaneous fit in each bin to determine the signal yield
 */

#ifndef DOUBLETAGYIELD
#define DOUBLETAGYIELD

#include"TTree.h"
#include"RooRealVar.h"
#include"RooFitResult.h"
#include"RooSimultaneous.h"
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
    /**
     * Plot projections of each bin in the fit
     */
    void PlotProjections(BinnedDataLoader *DataLoader, RooSimultaneous *Model);
    /**
     * Save signal yields from fit
     */
    void SaveSignalYields(const BinnedFitModel &FitModel, RooFitResult *Result, const Category &category) const;
  private:
    /**
     * The fit variable
     */
    RooRealVar m_SignalMBC;
    /**
     * Need this to apply a cut on the tag side
     */
    RooRealVar m_TagMBC;
    /**
     * The fit settings
     */
    Settings m_Settings;
    /**
     * TTree with double tag events
     */
    TTree *m_Tree;
};

#endif
