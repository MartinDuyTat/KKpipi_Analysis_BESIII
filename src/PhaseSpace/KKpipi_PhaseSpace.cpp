// Martin Duy Tat 26th November 2021

#include<algorithm>
#include<stdexcept>
#include"TTree.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"

KKpipi_PhaseSpace::KKpipi_PhaseSpace(TTree *Tree,
				     int Bins,
				     bool ReconstructedBins,
				     bool TrueBins): m_Momenta(16),
						     m_MomentaKalmanFit(16),
						     m_AmplitudePhaseSpace(Bins) {
  m_AmplitudePhaseSpace.SetBinEdges({1.20923});
  m_AmplitudePhaseSpace.UseVariableBinWidths(true);
  if(ReconstructedBins) {
    SetBranchAddresses_Rec(Tree);
  }
  if(TrueBins) {
    SetBranchAddresses_True(Tree);
  }
}

void KKpipi_PhaseSpace::SetBranchAddresses_Rec(TTree *Tree) {
  Tree->SetBranchAddress("SignalKalmanFitSuccess", &m_KalmanFitSuccess);
  Tree->SetBranchAddress("SignalKPluspx", m_Momenta.data() + 0);
  Tree->SetBranchAddress("SignalKPluspy", m_Momenta.data() + 1);
  Tree->SetBranchAddress("SignalKPluspz", m_Momenta.data() + 2);
  Tree->SetBranchAddress("SignalKPlusenergy", m_Momenta.data() + 3);
  Tree->SetBranchAddress("SignalKMinuspx", m_Momenta.data() + 4);
  Tree->SetBranchAddress("SignalKMinuspy", m_Momenta.data() + 5);
  Tree->SetBranchAddress("SignalKMinuspz", m_Momenta.data() + 6);
  Tree->SetBranchAddress("SignalKMinusenergy", m_Momenta.data() + 7);
  Tree->SetBranchAddress("SignalPiPluspx", m_Momenta.data() + 8);
  Tree->SetBranchAddress("SignalPiPluspy", m_Momenta.data() + 9);
  Tree->SetBranchAddress("SignalPiPluspz", m_Momenta.data() + 10);
  Tree->SetBranchAddress("SignalPiPlusenergy", m_Momenta.data() + 11);
  Tree->SetBranchAddress("SignalPiMinuspx", m_Momenta.data() + 12);
  Tree->SetBranchAddress("SignalPiMinuspy", m_Momenta.data() + 13);
  Tree->SetBranchAddress("SignalPiMinuspz", m_Momenta.data() + 14);
  Tree->SetBranchAddress("SignalPiMinusenergy", m_Momenta.data() + 15);
  Tree->SetBranchAddress("SignalKPluspxKalmanFit", m_MomentaKalmanFit.data() + 0);
  Tree->SetBranchAddress("SignalKPluspyKalmanFit", m_MomentaKalmanFit.data() + 1);
  Tree->SetBranchAddress("SignalKPluspzKalmanFit", m_MomentaKalmanFit.data() + 2);
  Tree->SetBranchAddress("SignalKPlusenergyKalmanFit", m_MomentaKalmanFit.data() + 3);
  Tree->SetBranchAddress("SignalKMinuspxKalmanFit", m_MomentaKalmanFit.data() + 4);
  Tree->SetBranchAddress("SignalKMinuspyKalmanFit", m_MomentaKalmanFit.data() + 5);
  Tree->SetBranchAddress("SignalKMinuspzKalmanFit", m_MomentaKalmanFit.data() + 6);
  Tree->SetBranchAddress("SignalKMinusenergyKalmanFit", m_MomentaKalmanFit.data() + 7);
  Tree->SetBranchAddress("SignalPiPluspxKalmanFit", m_MomentaKalmanFit.data() + 8);
  Tree->SetBranchAddress("SignalPiPluspyKalmanFit", m_MomentaKalmanFit.data() + 9);
  Tree->SetBranchAddress("SignalPiPluspzKalmanFit", m_MomentaKalmanFit.data() + 10);
  Tree->SetBranchAddress("SignalPiPlusenergyKalmanFit", m_MomentaKalmanFit.data() + 11);
  Tree->SetBranchAddress("SignalPiMinuspxKalmanFit", m_MomentaKalmanFit.data() + 12);
  Tree->SetBranchAddress("SignalPiMinuspyKalmanFit", m_MomentaKalmanFit.data() + 13);
  Tree->SetBranchAddress("SignalPiMinuspzKalmanFit", m_MomentaKalmanFit.data() + 14);
  Tree->SetBranchAddress("SignalPiMinusenergyKalmanFit", m_MomentaKalmanFit.data() + 15);
}

void KKpipi_PhaseSpace::SetBranchAddresses_True(TTree *Tree) {
  Tree->SetBranchAddress("NumberOfParticles", &m_TrueKinematics.NumberParticles);
  Tree->SetBranchAddress("ParticleIDs", m_TrueKinematics.ParticleIDs.data());
  Tree->SetBranchAddress("MotherIndex", m_TrueKinematics.MotherIndex.data());
  Tree->SetBranchAddress("True_Px", m_TrueKinematics.TruePx.data());
  Tree->SetBranchAddress("True_Py", m_TrueKinematics.TruePy.data());
  Tree->SetBranchAddress("True_Pz", m_TrueKinematics.TruePz.data());
  Tree->SetBranchAddress("True_Energy", m_TrueKinematics.TrueEnergy.data());
}
  
int KKpipi_PhaseSpace::KKpipiBin() const {
  return m_KalmanFitSuccess == 1 ? m_AmplitudePhaseSpace.WhichBin(m_MomentaKalmanFit) : m_AmplitudePhaseSpace.WhichBin(m_Momenta);
}

int KKpipi_PhaseSpace::TrueKKpipiBin() {
  std::vector<double> TrueMomenta;
  std::vector<int> DaughterIDs{321, -321, 211, -211};
  std::vector<int> IDs(m_TrueKinematics.ParticleIDs.begin(), m_TrueKinematics.ParticleIDs.begin() + m_TrueKinematics.NumberParticles);
  std::vector<double>::size_type SignalEnd_index;
  if(m_TrueKinematics.SignalD_index < m_TrueKinematics.TagD_index) {
    SignalEnd_index = m_TrueKinematics.TagD_index;
  } else {
    SignalEnd_index = m_TrueKinematics.NumberParticles;
  }
  for(auto DaughterID : DaughterIDs) {
    std::vector<double>::size_type Daughter_index = std::find(IDs.begin() + m_TrueKinematics.SignalD_index + 1, IDs.begin() + SignalEnd_index, DaughterID) - IDs.begin();
    TrueMomenta.push_back(m_TrueKinematics.TruePx[Daughter_index]);
    TrueMomenta.push_back(m_TrueKinematics.TruePy[Daughter_index]);
    TrueMomenta.push_back(m_TrueKinematics.TruePz[Daughter_index]);
    TrueMomenta.push_back(m_TrueKinematics.TrueEnergy[Daughter_index]);
  }
  return m_AmplitudePhaseSpace.WhichBin(TrueMomenta);
}

void KKpipi_PhaseSpace::FindDIndex() {
  // Copy the particle ID vector with the correct number of particles
  std::vector<int> IDs(m_TrueKinematics.ParticleIDs.begin(), m_TrueKinematics.ParticleIDs.begin() + m_TrueKinematics.NumberParticles);
  // Find the D0
  m_TrueKinematics.SignalD_index = std::find(IDs.begin(), IDs.end(), 421) - IDs.begin();
  // Find the D0bar
  m_TrueKinematics.TagD_index = std::find(IDs.begin(), IDs.end(), -421) - IDs.begin();
  // Find the D0 daughters and sort
  std::vector<int> D0Daughters(IDs.begin() + m_TrueKinematics.SignalD_index + 1, IDs.begin() + m_TrueKinematics.TagD_index);
  std::sort(D0Daughters.begin(), D0Daughters.end());
  // Find the D0bar daughters and sort
  std::vector<int> D0barDaughters(IDs.begin() + m_TrueKinematics.TagD_index + 1, IDs.end());
  std::sort(D0barDaughters.begin(), D0barDaughters.end());
  // State the KKpipi daughters and sort
  std::vector<int> SignalIDs{321, -321, 211, -211};
  std::sort(SignalIDs.begin(), SignalIDs.end());
  if(std::includes(D0Daughters.begin(), D0Daughters.end(), SignalIDs.begin(), SignalIDs.end())) {
    // If The D0 daughters include KKpipi, do nothing
  } else if(std::includes(D0barDaughters.begin(), D0barDaughters.end(), SignalIDs.begin(), SignalIDs.end())) {
    // If the D0bar daughters include KKpipi, swap the labels
    std::swap(m_TrueKinematics.SignalD_index, m_TrueKinematics.TagD_index);
  } else {
    throw std::logic_error("Could not find KKpipi truth information");
  }
}
