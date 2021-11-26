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
#include"DoubleTagYield.h"

class BinnedFitModel {
  public:
    /**
     * Constructor that takes in the signal MC TTree and sets up the signal shape for all the bins
     * @param settings Fit settings
     * @param Tree TTree with MC signal shape
     * @param SignalMBC The fit variable
     */
    BinnedFitModel(const Settings &settings, TTree *Tree, RooRealVar *SignalMBC);
    /**
     * Destructor that deletes the simultaneous PDF
     */
    ~BinnedFitModel();
    /**
     * Get the simultaneous PDF for double tag yield fit
     */
    RooSimultaneous* GetPDF();
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
    RooArgusBG *m_Argus = nullptr;
    /**
     * Initialize the signal component in a specific bin
     */
    void InitializeSignalShape(TTree *Tree);
    /**
     * Create the combinatorial component in a specific bin
     */
    void InitializeArgusShape();
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
    std::map<std::string, RooRealVar*> m_Yields;
};

#endif
