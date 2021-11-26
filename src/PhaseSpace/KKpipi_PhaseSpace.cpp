// Martin Duy Tat 26th November 2021

#include"TTree.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"

KKpipi_PhaseSpace::KKpipi_PhaseSpace(TTree *Tree, int Bins): m_Momenta(16), m_MomentaKalmanFit(16), m_AmplitudePhaseSpace(Bins) {
  m_AmplitudePhaseSpace.SetBinEdges({1.20923});
  m_AmplitudePhaseSpace.UseVariableBinWidths(true);
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

int KKpipi_PhaseSpace::KKpipiBin() const {
  return m_KalmanFitSuccess == 1 ? m_AmplitudePhaseSpace.WhichBin(m_MomentaKalmanFit) : m_AmplitudePhaseSpace.WhichBin(m_Momenta);
}
