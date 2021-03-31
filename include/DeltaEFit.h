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

class DeltaEFit {
  public:
    /**
     * Constructor that takes in a TTree and also initializes all parameters to standard values
     * @param Tree TTree object containing the data from BESIII
     * @param BranchName Name of the branch containing the \f$\Delta E\f$ variable
     */
    DeltaEFit(TTree *Tree);
    /**
     * Function for doing a fit of the \f$\Delta E\f$ distribution
     * First a binned fit with 200 bins is performed, then a more accurate unbinned fit is performed
     */
    void FitDeltaE();
    /**
     * Function for saving the fitted parameters to a text file
     * @param Filename Name of output text file
     */
    void SaveParameters(const std::string &Filename) const;
  private:
    /**
     * TTree containing the data we want to fit
     */
    TTree *m_Tree;
    /**
     * \f$\Delta E\f$ variable
     */
    RooRealVar m_DeltaE;
    /**
     * Signal yield of first Gaussian
     */
    RooRealVar m_Nsig1;
    /**
     * Signal yield of second Gaussian
     */
    RooRealVar m_Nsig2;
    /**
     * Combinatorial background yield
     */
    RooRealVar m_Nbkg;
    /**
     * Mean of first Gaussian
     */
    RooRealVar m_Mean1;
    /**
     * Mean of second Gaussian
     */
    RooRealVar m_Mean2;
    /**
     * Width of first Gaussian
     */
    RooRealVar m_Sigma1;
    /**
     * Width of second Gaussian
     */
    RooRealVar m_Sigma2;
    /**
     * Linear coefficient of polynomial
     */
    RooRealVar m_a;
    /**
     * Quadratic coefficient of polynomial
     */
    RooRealVar m_b;
};

#endif
