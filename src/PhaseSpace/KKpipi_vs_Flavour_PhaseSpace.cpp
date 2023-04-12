// Martin Duy Tat 26th November 2021

#include<utility>
#include<vector>
#include<stdexcept>
#include"TTree.h"
#include"PhaseSpace/KKpipi_vs_Flavour_PhaseSpace.h"

KKpipi_vs_Flavour_PhaseSpace::KKpipi_vs_Flavour_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, bool KSKK_binning): KKpipi_PhaseSpace(Tree, Bins, ReconstructedBins, TrueBins, KSKK_binning) {
  if(ReconstructedBins) {
    Tree->SetBranchAddress("TagKCharge", &m_KaonCharge);
  }
}

std::pair<int, int> KKpipi_vs_Flavour_PhaseSpace::Bin() const {
  return std::make_pair(-KKpipiBin()*m_KaonCharge, 0);
}

std::pair<int, int> KKpipi_vs_Flavour_PhaseSpace::TrueBin() {
  FindDIndex();
  std::vector<int>::size_type TagEnd_index;
  if(m_TrueKinematics.SignalD_index < m_TrueKinematics.TagD_index) {
    TagEnd_index = m_TrueKinematics.NumberParticles;
  } else {
    TagEnd_index = m_TrueKinematics.SignalD_index;
  }
  std::vector<int> IDs(m_TrueKinematics.ParticleIDs.begin(), m_TrueKinematics.ParticleIDs.begin() + m_TrueKinematics.NumberParticles);
  int KaonCharge;
  if(std::find(IDs.begin() + m_TrueKinematics.TagD_index + 1, IDs.begin() + TagEnd_index, 321) != IDs.begin() + TagEnd_index) {
    KaonCharge = +1;
  } else if(std::find(IDs.begin() + m_TrueKinematics.TagD_index + 1, IDs.begin() + TagEnd_index, -321) != IDs.begin() + TagEnd_index) {
    KaonCharge = -1;
  } else {
    throw std::logic_error("Cannot find true tag kaon ID");
  }
  return std::make_pair(-TrueKKpipiBin()*KaonCharge, 0);
}
