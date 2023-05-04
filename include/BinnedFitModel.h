// Martin Duy Tat 25th November 2021
/**
 * BinnedFitModel sets up a RooSimultaneous model for the binned fit of a given tag mode
 */

#ifndef BINNEDFITMODEL
#define BINNEDFITMODEL

#include<string>
#include<map>
#include"TTree.h"
#include"RooFFTConvPdf.h"
#include"RooSimultaneous.h"
#include"RooArgusBG.h"
#include"RooAddPdf.h"
#include"RooAbsPdf.h"
#include"Settings.h"
#include"Category.h"
#include"CholeskySmearing.h"
#include"RooShapes/FitShape.h"

class DoubleTagYield;

class BinnedFitModel {
  public:
    /**
     * Constructor that takes in the signal MC TTree and sets up the signal shape for all the bins
     * @param settings Fit settings
     * @param SignalMBC The fit variable
     */
    BinnedFitModel(const Settings &settings, RooRealVar *SignalMBC);
    /**
     * Destructor that deletes the simultaneous PDF and the peaking background shapes
     */
    ~BinnedFitModel();
    /**
     * Get the simultaneous PDF for double tag yield fit
     */
    RooSimultaneous* GetPDF();
    /**
     * Reset the signal yields to generator values
     */
    void SetGeneratorYields();
    /**
     * DoubleTagYield is a friend so that it can access the yield variables
     */
    friend class DoubleTagYield;
  private:
    /**
     * The fit variable
     */
    RooRealVar *m_SignalMBC;
    /**
     * Simultaneous PDF for the binned fit
     */
    RooSimultaneous *m_PDF = nullptr;
    /**
     * The signal shape from signal MC convolved with a double Gaussian
     */
    RooFFTConvPdf *m_SignalShapeConv = nullptr;
    /**
     * The combinatorial shape which floats in the fit, but shared between all bins
     */
    RooAbsPdf *m_Combinatorial = nullptr;
    /**
     * Map of all the peaking background shapes
     */
    std::map<std::string, FitShape*> m_PeakingBackgroundShapes;
    /**
     * Initialize all the yield variables
     */
    void InitializeYields();
    /**
     * Initialize the signal component for all bins
     */
    void InitializeSignalShape();
    /**
     * Create the combinatorial component for all bins
     */
    void InitializeCombinatorialShape();
    /**
     * Create the peaking background shape components for all bins
     */
    void InitializePeakingBackgroundShapes();
    /**
     * Create the total PDF and yield in a specific bin
     */
    RooAddPdf* CreateBinPDF(const std::string &CategoryString);
    /**
     * Initialize the full PDF in all bins and add them to the simultaneous fit
     */
    void InitializePDF();
    /**
     * The object keeping track of all the categories
     */
    Category m_Category;
    /**
     * The fit settings
     */
    Settings m_Settings;
    /**
     * Map of all yield parameters
     */
    std::map<std::string, RooAbsReal*> m_Yields;
    /**
     * Map of Gaussian resolution and Argus parameters
     */
    std::map<std::string, RooRealVar*> m_Parameters;
    /**
     * The signal yields that MINOS will run over
     */
    RooArgSet m_SignalYields;
    /**
     * Smear the peaking backgrounds to estimate the systematic uncertainties
     */
    void SmearPeakingBackgrounds();
    /**
     * Map of the Cholesky smearing objects
     */
    std::map<std::string, CholeskySmearing> m_CholeskyDecompositions;
    /**
     * Helper function that sets up all the Cholesky decompositions for smearing of correlated peaking backgrounds
     */
    void PrepareSmearing();
};

#endif
