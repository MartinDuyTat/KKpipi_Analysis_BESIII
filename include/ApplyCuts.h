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
     * @param TTree or TChain with event
     */
    TTree* operator()(TTree *InputTree) const;
  private:
    /**
     * Cuts that will be applied to selection
     */
    TCut m_Cuts;
};

#endif
