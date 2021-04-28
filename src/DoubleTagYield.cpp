// Martin Duy Tat 28th April 2021

#include<vector>
#include"TTree.h"
#include"TCut.h"
#include"TEntryList.h"
#include"TMath.h"
#include"AmplitudePhaseSpace.h"
#include"Event.h"
#include"DoubleTagYield.h"

DoubleTagYield::DoubleTagYield(int NBins, TCut Cuts): m_NBins(NBins), m_AmplitudePhaseSpace(AmplitudePhaseSpace(m_NBins)), m_BinYieldS(std::vector<int>(m_NBins)), m_BinYieldA(std::vector<int>(m_NBins)), m_BinYieldB(std::vector<int>(m_NBins)), m_BinYieldC(std::vector<int>(m_NBins)), m_BinYieldD(std::vector<int>(m_NBins)), m_EventsOutsidePhaseSpace(0), m_EventsOutsideMBCSpace(0), m_Cuts(Cuts) {
}

char DoubleTagYield::DetermineMBCRegion(double SignalMBC, double TagMBC) const {
  if(SignalMBC > 1.86 && SignalMBC < 1.87 && TagMBC > 1.86 && TagMBC < 1.87) {
    return 'S';
  } else if(SignalMBC > 1.86 && SignalMBC < 1.87 && TagMBC > 1.83 && TagMBC < 1.855) {
    return 'A';
  } else if(SignalMBC > 1.83 && SignalMBC < 1.855 && TagMBC > 1.86 && TagMBC < 1.87) {
    return 'B';
  } else if(SignalMBC > 1.83 && SignalMBC < 1.855 && TagMBC > 1.83 && TagMBC < 1.855) {
    if(TMath::Abs(SignalMBC - TagMBC) < 0.0035) {
      return 'C';
    } else if(TMath::Abs(SignalMBC - TagMBC) > 0.0055) {
      return 'D';
    } else {
      return 'F';
    }
  } else {
    return 'F';
  }
}

void DoubleTagYield::SetDaughterBranchAddresses(TTree *Tree, std::vector<double> &Momenta, std::vector<double> &MomentaKalmanFit) const {
  Momenta = std::vector<double>(16);
  MomentaKalmanFit = std::vector<double>(16);
  Tree->SetBranchAddress("SignalKPluspx", Momenta.data() + 0);
  Tree->SetBranchAddress("SignalKPluspy", Momenta.data() + 1);
  Tree->SetBranchAddress("SignalKPluspz", Momenta.data() + 2);
  Tree->SetBranchAddress("SignalKPlusenergy", Momenta.data() + 3);
  Tree->SetBranchAddress("SignalKMinuspx", Momenta.data() + 4);
  Tree->SetBranchAddress("SignalKMinuspy", Momenta.data() + 5);
  Tree->SetBranchAddress("SignalKMinuspz", Momenta.data() + 6);
  Tree->SetBranchAddress("SignalKMinusenergy", Momenta.data() + 7);
  Tree->SetBranchAddress("SignalPiPluspx", Momenta.data() + 8);
  Tree->SetBranchAddress("SignalPiPluspy", Momenta.data() + 9);
  Tree->SetBranchAddress("SignalPiPluspz", Momenta.data() + 10);
  Tree->SetBranchAddress("SignalPiPlusenergy", Momenta.data() + 11);
  Tree->SetBranchAddress("SignalPiMinuspx", Momenta.data() + 12);
  Tree->SetBranchAddress("SignalPiMinuspy", Momenta.data() + 13);
  Tree->SetBranchAddress("SignalPiMinuspz", Momenta.data() + 14);
  Tree->SetBranchAddress("SignalPiMinusenergy", Momenta.data() + 15);
  Tree->SetBranchAddress("SignalKPluspxKalmanFit", MomentaKalmanFit.data() + 0);
  Tree->SetBranchAddress("SignalKPluspyKalmanFit", MomentaKalmanFit.data() + 1);
  Tree->SetBranchAddress("SignalKPluspzKalmanFit", MomentaKalmanFit.data() + 2);
  Tree->SetBranchAddress("SignalKPlusenergyKalmanFit", MomentaKalmanFit.data() + 3);
  Tree->SetBranchAddress("SignalKMinuspxKalmanFit", MomentaKalmanFit.data() + 4);
  Tree->SetBranchAddress("SignalKMinuspyKalmanFit", MomentaKalmanFit.data() + 5);
  Tree->SetBranchAddress("SignalKMinuspzKalmanFit", MomentaKalmanFit.data() + 6);
  Tree->SetBranchAddress("SignalKMinusenergyKalmanFit", MomentaKalmanFit.data() + 7);
  Tree->SetBranchAddress("SignalPiPluspxKalmanFit", MomentaKalmanFit.data() + 8);
  Tree->SetBranchAddress("SignalPiPluspyKalmanFit", MomentaKalmanFit.data() + 9);
  Tree->SetBranchAddress("SignalPiPluspzKalmanFit", MomentaKalmanFit.data() + 10);
  Tree->SetBranchAddress("SignalPiPlusenergyKalmanFit", MomentaKalmanFit.data() + 11);
  Tree->SetBranchAddress("SignalPiMinuspxKalmanFit", MomentaKalmanFit.data() + 12);
  Tree->SetBranchAddress("SignalPiMinuspyKalmanFit", MomentaKalmanFit.data() + 13);
  Tree->SetBranchAddress("SignalPiMinuspzKalmanFit", MomentaKalmanFit.data() + 14);
  Tree->SetBranchAddress("SignalPiMinusenergyKalmanFit", MomentaKalmanFit.data() + 15);
}

void DoubleTagYield::CalculateBinnedRawYields(TTree *Tree) {
  Tree->Draw(">> elist", m_Cuts, "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  Tree->SetEntryList(elist);
  double SignalMBC, TagMBC;
  int KalmanFitSuccess;
  Tree->SetBranchAddress("SignalMBC", &SignalMBC);
  Tree->SetBranchAddress("TagMBC", &TagMBC);
  Tree->SetBranchAddress("SignalKalmanFitSuccess", &KalmanFitSuccess);
  std::vector<double> Momenta, MomentaKalmanFit;
  SetDaughterBranchAddresses(Tree, Momenta, MomentaKalmanFit);
  for(int i = 0; i < elist->GetN(); i++) {
    Tree->GetEntry(Tree->GetEntryNumber(i));
    char Region = DetermineMBCRegion(SignalMBC, TagMBC);
    if(Region == 'F') {
      m_EventsOutsideMBCSpace++;
      continue;
    }
    int BinNumber = KalmanFitSuccess == 1 ? m_AmplitudePhaseSpace.WhichBin(Event(MomentaKalmanFit)) : m_AmplitudePhaseSpace.WhichBin(Event(Momenta));
    int BinIndex = TMath::Abs(BinNumber) - 1;
    if(BinNumber == 0) {
      m_EventsOutsidePhaseSpace++;
      continue;
    }
    if(Region == 'S') {
      m_BinYieldS[BinIndex]++;
    } else if(Region == 'A') {
      m_BinYieldS[BinIndex]++;
    } else if(Region == 'B') {
      m_BinYieldB[BinIndex]++;
    } else if(Region == 'C') {
      m_BinYieldC[BinIndex]++;
    } else if(Region == 'D') {
      m_BinYieldD[BinIndex]++;
    }
  }
}

double DoubleTagYield::GetBinYield(int i) const {
  int BinIndex = TMath::Abs(i) - 1;
  // Area in beam constrained 2D plane
  double A_S = 10*10;
  double A_A = 10*25;
  double A_B = 25*10;
  double A_D = 19.5*19.5;
  double A_C = 25*25 - 21.5*21.5;
  double Background = (A_S/A_D)*m_BinYieldD[BinIndex] + (A_S/A_A)*(m_BinYieldA[BinIndex] - (A_S/A_A)*m_BinYieldD[BinIndex]) + (A_S/A_B)*(m_BinYieldB[BinIndex] - (A_S/A_B)*m_BinYieldD[BinIndex]) + (A_S/A_C)*(m_BinYieldC[BinIndex] - (A_S/A_C)*m_BinYieldD[BinIndex]);
  return m_BinYieldS[BinIndex] - Background;
}

int DoubleTagYield::GetEventsOutsidePhaseSpace() const {
  return m_EventsOutsidePhaseSpace;
}

int DoubleTagYield::GetEventsOutsideMBCSpace() const {
  return m_EventsOutsideMBCSpace;
}
