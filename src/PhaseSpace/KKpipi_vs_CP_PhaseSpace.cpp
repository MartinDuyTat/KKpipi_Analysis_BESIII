// Martin Duy Tat 21st December 2021

#include<utility>
#include"TTree.h"
#include"TMath.h"
#include"PhaseSpace/KKpipi_vs_CP_PhaseSpace.h"

KKpipi_vs_CP_PhaseSpace::KKpipi_vs_CP_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, bool KSKK_binning): KKpipi_PhaseSpace(Tree, Bins, ReconstructedBins, TrueBins, KSKK_binning) {
}

std::pair<int, int> KKpipi_vs_CP_PhaseSpace::Bin() const {
  return std::make_pair(TMath::Abs(KKpipiBin()), 0);
}

std::pair<int, int> KKpipi_vs_CP_PhaseSpace::TrueBin() {
  FindDIndex();
  return std::make_pair(TMath::Abs(TrueKKpipiBin()), 0);
}
