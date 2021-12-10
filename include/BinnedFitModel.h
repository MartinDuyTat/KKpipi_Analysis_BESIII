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
#include"Settings.h"
#include"Category.h"

class DoubleTagYield;

class BinnedFitModel {
  public:
    /**
     * Constructor that takes in the signal MC TTree and sets up the signal shape for all the bins
     * @param settings Fit settings
     * @param SignalMBC The fit variable
     */
    BinnedFitModel(const Settings &settings, RooRealVar *SignalMBC, RooRealVar *TagMBC);
    /**
     * Destructor that deletes the simultaneous PDF
     */
    ~BinnedFitModel();
    /**
     * Get the simultaneous PDF for double tag yield fit
     */
    RooSimultaneous* GetPDF();
    /**
     * Get the fraction of events inside the signal region
     */
    double GetFractionInSignalRegion() const;
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
     * Need this to apply a cut on the tag side
     */
    RooRealVar *m_TagMBC;
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
    RooArgusBG *m_Argus = nullptr;
    /**
     * Map of all the peaking background shapes
     */
    std::map<std::string, RooAbsPdf*> m_PeakingBackgroundShapes;
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
    void InitializeArgusShape();
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
};

#endif
