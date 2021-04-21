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
#include"TTree.h"
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooExtendPdf.h"

class SingleTagYield {
  public:
    /**
     * Constructor that takes in the TTree object with the data and a TTree with the MC signal shape
     * @param DataTree TTree with signal events
     * @param MCSignalTree TTree with MC signal shape for fitting the signal
     * @param TreeName Name of the TTree
     */
    SingleTagYield(TTree *DataTree, TTree *MCSignalTree);
    /**
     * Destructor that deletes all the PDFs and corresponding variables for peaking backgrounds that were allocated on the heap
     */
    ~SingleTagYield();
    /**
     * Function that parses peaking background components from a file and places a fixed component for this in the fit
     * Text file must have each peaking background on a separate line, in the format "Name Mean Sigma Yield"
     * @param Filename Filename of text file with peaking background Gaussian components
     */
    void AddPeakingComponent(const std::string &Filename);
    /**
     * Function that performs the fit of the single tag yield in RooFit
     * @param TagMode Plot label tag mode (#pi#pi for pipi etc)
     * @param Filename Filename to save plot to
     */
    void FitYield(const std::string &TagMode, const std::string &Filename);
    /**
     * Function that saves the fit parameters to a text file
     * @param Filename Filename of text file
     */
    void SaveFitParameters(const std::string &Filename) const;
    /**
     * Function that returns fitted yield, only call this after fit
     */
    double GetYield() const;
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
     * Beam constrained mass variable, independent parameter in fit
     */
    RooRealVar m_MBC;
    /**
     * An Argus PDF shape parameter
     */
    RooRealVar m_c;
    /**
     * Endpoint of Argus PDF shape
     */
    RooRealVar m_End;
    /**
     * Mean of Gaussian modelling detector resolution
     */
    RooRealVar m_Mean;
    /**
     * Width of Gaussian modelling detector resolution
     */
    RooRealVar m_Sigma;
    /**
     * Signal yield
     */
    RooRealVar m_Nsig;
    /**
     * Background yield
     */
    RooRealVar m_Nbkg;
    /**
     * Weighting to account for luminosity scale
     */
    RooRealVar m_LuminosityWeight;
    /**
     * Vector of names of peaking background
     */
    std::vector<std::string> m_PeakingName;
    /**
     * Vector of RooRealVar objects for mean of peaking background Gaussians
     */
    std::vector<RooRealVar*> m_PeakingMean;
    /**
     * Vector of RooRealVar objects for width of peaking background Gaussians
     */
    std::vector<RooRealVar*> m_PeakingSigma;
    /**
     * Vector of RooRealVar objects for yield of peaking background Gaussians
     */
    std::vector<RooRealVar*> m_PeakingYield;
    /**
     * Vector of RooGaussian objects for PDF of peaking backgrounds
     */
    std::vector<RooGaussian*> m_PeakingPDF;
    /**
     * Vector of RooExtendPDF objects for PDF of peaking backgrounds
     */
    std::vector<RooExtendPdf*> m_PeakingExPDF;
};

#endif
