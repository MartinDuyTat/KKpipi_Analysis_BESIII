// Martin Duy Tat 25th February 2021
/**
 * ApplyCuts is a class that applies cuts to a TTree or TChain and outputs a TTree with the events that passed the selection
 * It is used as a functor after the cuts are specified in the constructor
 */

#ifndef APPLYCUTS
#define APPLYCUTS

#include"TCut.h"
#include"TTree.h"

class ApplyCuts {
  public:
    /**
     * Constructor that takes in a cut for the selection
     * @param Cuts Cuts that will be applied in selection
     */
    ApplyCuts(const TCut &Cuts);
    /**
     * () operator overload so that one can pass a TTree or TChain to this object and get a skimmed TTree back
     * If a dataset type and luminosity scale is given, these branches are also added to the final TTree
     * @param TTree or TChain with event
     * @param DataSetType Integer between \f$0\f$ and \f$9\f$, labelling the dataset (description in PrepareTagTree application)
     * @param LuminosityScale Luminosity scale of MC, TTree will be filled with the inverse of this to scale MC to that of data
     */
    TTree* operator()(TTree *InputTree, int DataSetType = -1, double LuminosityScale = 1.0) const;
  private:
    /**
     * Cuts that will be applied to selection
     */
    TCut m_Cuts;
};

#endif
