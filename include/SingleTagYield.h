// Martin Duy Tat 4th April 2021
/** 
 * SingleTagYield is a class for fitting the single tag yield
 * The yield is fitted with a signal shape from an exclusive signal MC sample and an Argus PDF for the combinatorial background
 * Any peaking backgrounds are fixed by Gaussians
 */

#ifndef SINGLETAGYIELD
#define SINGLETAGYIELD

#include<string>
#include"TTree.h"
#include"RooRealVar.h"

class SingleTagYield {
  public:
    /**
     * Constructor that takes in the TTree object with the data and a TTree with the MC signal shape
     * @param DataTree TTree with signal events
     * @param MCSignalTree TTree with MC signal shape for fitting the signal
     * @param TreeName Name of the TTree
     */
    SingleTagYield(TTree *DataTree, TTree *MCSignalTree);
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
     * Signal yield
     */
    RooRealVar m_Nsig;
    /**
     * Background yield
     */
    RooRealVar m_Nbkg;
};

#endif
