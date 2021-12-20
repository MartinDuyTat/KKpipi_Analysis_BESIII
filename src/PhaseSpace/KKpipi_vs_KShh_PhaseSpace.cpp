// Martin Duy Tat 20th December 2021

#include<string>
#include<stdexcept>
#include<utility>
#include<vector>
#include<algorithm>
#include"TFile.h"
#include"TTree.h"
#include"TLorentzVector.h"
#include"PhaseSpace/KKpipi_vs_KShh_PhaseSpace.h"

KKpipi_vs_KShh_PhaseSpace::KKpipi_vs_KShh_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, const std::string &Mode): KKpipi_PhaseSpace(Tree, Bins, ReconstructedBins, TrueBins), m_Mode(Mode) {
  std::string BinningFilename;
  if(Mode == "KSpipi") {
    BinningFilename = std::string(BINNING_SCHEME_DIR) + "KsPiPi_equal.root";
    Tree->SetBranchAddress("TagPiPluspx", &m_hPlus_P[0]);
    Tree->SetBranchAddress("TagPiPluspy", &m_hPlus_P[1]);
    Tree->SetBranchAddress("TagPiPluspz", &m_hPlus_P[2]);
    Tree->SetBranchAddress("TagPiPlusenergy", &m_hPlus_P[3]);
    Tree->SetBranchAddress("TagPiPluspxKalmanFit", &m_hPlus_P_KalmanFit[0]);
    Tree->SetBranchAddress("TagPiPluspyKalmanFit", &m_hPlus_P_KalmanFit[1]);
    Tree->SetBranchAddress("TagPiPluspzKalmanFit", &m_hPlus_P_KalmanFit[2]);
    Tree->SetBranchAddress("TagPiPlusenergyKalmanFit", &m_hPlus_P_KalmanFit[3]);
    Tree->SetBranchAddress("TagPiMinuspx", &m_hMinus_P[0]);
    Tree->SetBranchAddress("TagPiMinuspy", &m_hMinus_P[1]);
    Tree->SetBranchAddress("TagPiMinuspz", &m_hMinus_P[2]);
    Tree->SetBranchAddress("TagPiMinusenergy", &m_hMinus_P[3]);
    Tree->SetBranchAddress("TagPiMinuspxKalmanFit", &m_hMinus_P_KalmanFit[0]);
    Tree->SetBranchAddress("TagPiMinuspyKalmanFit", &m_hMinus_P_KalmanFit[1]);
    Tree->SetBranchAddress("TagPiMinuspzKalmanFit", &m_hMinus_P_KalmanFit[2]);
    Tree->SetBranchAddress("TagPiMinusenergyKalmanFit", &m_hMinus_P_KalmanFit[3]);
  /*} else if(Mode == "KSKK") {
    BinningFilename = std::string(BINNING_SCHEME_DIR) + "KsKK_2bins.root";*/
  } else {
    throw std::invalid_argument("Mode " + Mode + " is not a valid KShh mode");
  }
  TFile BinningFile(BinningFilename.c_str(), "READ");
  BinningFile.GetObject("dkpp_bin_h", m_BinningScheme);
  m_BinningScheme->SetDirectory(0);
  BinningFile.Close();
  Tree->SetBranchAddress("TagKalmanFitSuccess", &m_KalmanFitSuccess);
  Tree->SetBranchAddress("TagKSpx", &m_KS_P[0]);
  Tree->SetBranchAddress("TagKSpy", &m_KS_P[1]);
  //Tree->SetBranchAddress("TagKSpz", &m_KS_P[2]);
  Tree->SetBranchAddress("TagKSz", &m_KS_P[2]); // Typo in current BOSS version, need to rerun!
  Tree->SetBranchAddress("TagKSenergy", &m_KS_P[3]);
  Tree->SetBranchAddress("TagKSpxKalmanFit", &m_KS_P_KalmanFit[0]);
  Tree->SetBranchAddress("TagKSpyKalmanFit", &m_KS_P_KalmanFit[1]);
  Tree->SetBranchAddress("TagKSpzKalmanFit", &m_KS_P_KalmanFit[2]);
  Tree->SetBranchAddress("TagKSenergyKalmanFit", &m_KS_P_KalmanFit[3]);
}

std::pair<int, int> KKpipi_vs_KShh_PhaseSpace::Bin() const {
  return std::make_pair(KKpipiBin(), GetKShhBin());
}

std::pair<int, int> KKpipi_vs_KShh_PhaseSpace::TrueBin() {
  return std::make_pair(TrueKKpipiBin(), GetTrueKShhBin());
}

int KKpipi_vs_KShh_PhaseSpace::LookUpBinNumber(double M2Plus, double M2Minus) const {
  Float_t BinNumberFloat = m_BinningScheme->GetBinContent(m_BinningScheme->GetXaxis()->FindBin(M2Plus), m_BinningScheme->GetYaxis()->FindBin(M2Minus));
  return static_cast<int>(BinNumberFloat);
}

int KKpipi_vs_KShh_PhaseSpace::GetMappedKShhBin(double M2Plus, double M2Minus) const {
  int x = m_BinningScheme->GetXaxis()->FindBin(M2Plus);
  int y = m_BinningScheme->GetYaxis()->FindBin(M2Minus);
  int i = 1;
  // Loop in circles around this point until reaching the Dalitz boundary
  while(true) {
    if(i > 1000) {
      throw std::runtime_error("Dalitz point too far outside phase space");
    }
    for(int j = 0; j <= i; j++) {
      // All possible combinations of displacements at the same distance
      std::vector<std::pair<int, int>> BinList{{i, j}, {i, -j}, {-i, j}, {-i, -j}, {j, i}, {j, -i}, {-j, i}, {-j, -i}};
      for(auto iter = BinList.begin(); iter != BinList.end(); iter++) {
        int NewBin = static_cast<int>(m_BinningScheme->GetBinContent(x + iter->first, y + iter->second));
        // Once we reach the Dalitz boundary the bin number is non-zero
        if(NewBin != 0) {
          return NewBin;
        }
      }
    }
    i++;
  }
}

int KKpipi_vs_KShh_PhaseSpace::GetKShhBin() const {
  double M2Plus = m_KalmanFitSuccess == 1 ? (m_KS_P_KalmanFit + m_hPlus_P_KalmanFit).M2() : (m_KS_P + m_hPlus_P).M2();
  double M2Minus = m_KalmanFitSuccess == 1 ? (m_KS_P_KalmanFit + m_hMinus_P_KalmanFit).M2() : (m_KS_P + m_hMinus_P).M2();
  int Bin = LookUpBinNumber(M2Plus, M2Minus);
  if(Bin != 0) {
    return Bin;    
  } else {
    return GetMappedKShhBin(M2Plus, M2Minus);

  }
}

int KKpipi_vs_KShh_PhaseSpace::GetTrueKShhBin() {
  // Get index of signal and tag D in truth information
  FindDIndex();
  // Get the end index of the tag D
  std::vector<int>::size_type TagEnd_index;
  if(m_TrueKinematics.SignalD_index < m_TrueKinematics.TagD_index) {
    TagEnd_index = m_TrueKinematics.NumberParticles;
  } else {
    TagEnd_index = m_TrueKinematics.SignalD_index;
  }
  // Indices of KS, hPlus, hMinus in truth information
  std::vector<int>::size_type KS_index, hPlus_index, hMinus_index;
  // Particle ID of h+ and h-
  int h_ID = m_Mode == "KSpipi" ? 211 : 321;
  // begin() iterator of true IDs
  auto ID_begin = m_TrueKinematics.ParticleIDs.begin();
  // Get index of KS in truth information
  auto KS_iter = std::find(ID_begin + m_TrueKinematics.TagD_index + 1, ID_begin + TagEnd_index, 310);
  if(KS_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find KS in truth information");
  } else {
    KS_index = KS_iter - ID_begin;
  }
  // Repeat for hPlus
  auto hPlus_iter = std::find(KS_iter + 3, ID_begin + TagEnd_index, h_ID);
  if(hPlus_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find h+ in truth information");
  } else {
    hPlus_index = hPlus_iter - ID_begin;
  }
  // Repeat for hMinus
  auto hMinus_iter = std::find(KS_iter + 3, ID_begin + TagEnd_index, -h_ID);
  if(hMinus_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find h- in truth information");
  } else {
    hMinus_index = hMinus_iter - ID_begin;
  }
  // Calculate Dalitz coordinates
  TLorentzVector KS_P(m_TrueKinematics.TruePx[KS_index], m_TrueKinematics.TruePy[KS_index], m_TrueKinematics.TruePz[KS_index], m_TrueKinematics.TrueEnergy[KS_index]);
  TLorentzVector hPlus_P(m_TrueKinematics.TruePx[hPlus_index], m_TrueKinematics.TruePy[hPlus_index], m_TrueKinematics.TruePz[hPlus_index], m_TrueKinematics.TrueEnergy[hPlus_index]);
  TLorentzVector hMinus_P(m_TrueKinematics.TruePx[hMinus_index], m_TrueKinematics.TruePy[hMinus_index], m_TrueKinematics.TruePz[hMinus_index], m_TrueKinematics.TrueEnergy[hMinus_index]);
  double M2Plus = (KS_P + hPlus_P).M2();
  double M2Minus = (KS_P + hMinus_P).M2();
  // Get bin number
  int Bin = LookUpBinNumber(M2Plus, M2Minus);
  if(Bin != 0) {
    return Bin;
  } else {
    return GetMappedKShhBin(M2Plus, M2Minus);
  }
}
