// Martin Duy Tat 31st March 2021
/**
 * DeltaEFit is a class for doing a fit on \f$\Delta E\f$ in an input TTree
 * A double Gaussian with a polynomial background is fitted
 */

#ifndef DELTAEFIT
#define DELTAEFIT

#include<string>
#include"TTree.h"
#include"RooRealVar.h"
#include"RooDataHist.h"
#include"RooDataSet.h"
#include"RooAbsPdf.h"
#include"RooFitResult.h"
#include"Settings.h"
#include"DeltaEFitModel.h"

class DeltaEFit {
  public:
    /**
     * Constructor that takes in a TTree and also initializes all parameters to standard values
     * @param Tree TTree object containing the data from BESIII
     * @param settings The fit settings
     */
    DeltaEFit(TTree *Tree, const Settings &settings);
    /**
     * Function for doing a fit of the \f$\Delta E\f$ distribution
     * First a binned fit with 1000 bins is performed, then a more accurate unbinned fit is performed if necessary
     */
    void FitDeltaE();
    /**
     * Save the plot
     * @param UnbinnedData The dataset
     * @param Model The fit model after fit
     */
    void SavePlot(const RooDataSet &UnbinnedData, RooAbsPdf *Model) const;
    /**
     * Function for saving the fitted parameters to a text file
     * @param Results The fit results
     */
    void SaveParameters(RooFitResult *Results);
  private:
    /**
     * TTree containing the data we want to fit
     */
    TTree *m_Tree;
    /**
     * Fit settings
     */
    Settings m_Settings;
    /**
     * \f$\Delta E\f$ variable
     */
    RooRealVar m_DeltaE;
    /**
     * Weighting to account for luminosity scale
     */
    RooRealVar m_LuminosityWeight;
    /**
     * The lower \f$\Delta E\f$ cut
     */
    double m_DeltaE_Low;
    /**
     * The upper \f$\Delta E\f$ cut
     */
    double m_DeltaE_High;
};

#endif
