// Martin Duy Tat 4th April 2021

#include<string>
#include"TTree.h"
#include"RooRealVar.h"
#include"SingleTagYield.h"

SingleTagYield::SingleTagYield(TTree *DataTree, TTree *MCSignalTree):
                               m_DataTree(DataTree),
			       m_MCSignalTree(MCSignalTree),
			       m_MBC(RooRealVar("MBC", "MBC", 1.83, 1.89)),
			       m_c(RooRealVar("c", "c", -10, -100, 100)),
			       m_End(RooRealVar("End", "End", 1.8865)),
			       m_Nsig(RooRealVar("Nsig", "Nsig", m_DataTree->GetEntries()/2, 0.0, m_DataTree->GetEntries())),
			       m_Nbkg(RooRealVar("Nbkg", "Nbkg", m_DataTree->GetEntries()/2, 0.0, m_DataTree->GetEntries())) {
}
