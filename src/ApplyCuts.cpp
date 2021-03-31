// Martin Duy Tat 25th February 2021

#include"ApplyCuts.h"
#include"TCut.h"
#include"TTree.h"
#include"TEntryList.h"
#include"TDirectory.h"

ApplyCuts::ApplyCuts(const TCut &Cuts): m_Cuts(Cuts) {
}

TTree* ApplyCuts::operator()(TTree *InputTree) const {
  InputTree->Draw(">> elist", m_Cuts, "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  InputTree->SetEntryList(elist);
  TTree *OutputTree = InputTree->CloneTree(0);
  for(Long64_t i = 0; i < elist->GetN(); i++) {
    InputTree->GetEntry(InputTree->GetEntryNumber(i));
    OutputTree->Fill();
  }
  return OutputTree;
}
