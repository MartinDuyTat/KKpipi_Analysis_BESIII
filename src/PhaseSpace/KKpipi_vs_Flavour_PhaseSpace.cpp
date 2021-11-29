// Martin Duy Tat 26th November 2021

#include<utility>
#include"TTree.h"
#include"PhaseSpace/KKpipi_vs_Flavour_PhaseSpace.h"

KKpipi_vs_Flavour_PhaseSpace::KKpipi_vs_Flavour_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins): KKpipi_PhaseSpace(Tree, Bins, ReconstructedBins, TrueBins) {
  Tree->SetBranchAddress("TagKCharge", &m_KaonCharge);
}

std::pair<int, int> KKpipi_vs_Flavour_PhaseSpace::Bin() const {
  return std::make_pair(KKpipiBin()*m_KaonCharge, 0);
}

std::pair<int, int> KKpipi_vs_Flavour_PhaseSpace::TrueBin() const {
  return std::make_pair(TrueKKpipiBin(), 0);
}
