// Martin Duy Tat 20th December 2021

#include<string>
#include<stdexcept>
#include<utility>
#include<vector>
#include<algorithm>
#include"TFile.h"
#include"TTree.h"
#include"TLorentzVector.h"
#include"PhaseSpace/KKpipi_vs_K0hh_PhaseSpace.h"
#include"PhaseSpace/DalitzUtilities.h"

KKpipi_vs_K0hh_PhaseSpace::KKpipi_vs_K0hh_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, const std::string &Mode, bool KSKK_binning, bool KKpipiPartReco, bool KStoKLBackground): KKpipi_PhaseSpace(Tree, Bins, ReconstructedBins, TrueBins, KSKK_binning), m_K0Mode(Mode.substr(0, 2)), m_hhMode(Mode.substr(2)), m_KStoKLBackground(KStoKLBackground) {
  if(m_K0Mode != "KS" && m_K0Mode != "KL") {
    throw std::invalid_argument("K0 mode " + m_K0Mode + " not valid");
  } 
  if(m_hhMode != "pipi" && m_hhMode != "KK") {
    throw std::invalid_argument("hh mode " + m_hhMode + " not valid");
  }
  std::string BinningFilename;
  if(m_hhMode == "pipi") {
    BinningFilename = std::string(BINNING_SCHEME_DIR) + "KsPiPi_equal.root";
    if(ReconstructedBins) {
      SetK0hhBranchAddresses(Tree);
    }
  } else {
    BinningFilename = std::string(BINNING_SCHEME_DIR) + "KsKK_2bins.root";
  }
  TFile BinningFile(BinningFilename.c_str(), "READ");
  BinningFile.GetObject("dkpp_bin_h", m_BinningScheme);
  m_BinningScheme->SetDirectory(0);
  BinningFile.Close();
  if(ReconstructedBins) {
    SetK0hhBranchAddresses(Tree);
    if(!KKpipiPartReco) {
      Tree->SetBranchAddress("TagKalmanFitSuccess", &m_KalmanFitSuccess);
    } else {
      if(m_K0Mode == "KL") {
	throw std::invalid_argument("Cannot have two missing particles");
      }
      Tree->SetBranchAddress("SignalKalmanFitSuccess", &m_KalmanFitSuccess);
    }
  }
}

void KKpipi_vs_K0hh_PhaseSpace::SetK0hhBranchAddresses(TTree *Tree) {
  const std::string hName = m_hhMode == "pipi" ? "Pi" : "K";
  Tree->SetBranchAddress(("Tag" + hName + "Pluspx").c_str(), &m_hPlus_P[0]);
  Tree->SetBranchAddress(("Tag" + hName + "Pluspy").c_str(), &m_hPlus_P[1]);
  Tree->SetBranchAddress(("Tag" + hName + "Pluspz").c_str(), &m_hPlus_P[2]);
  Tree->SetBranchAddress(("Tag" + hName + "Plusenergy").c_str(), &m_hPlus_P[3]);
  Tree->SetBranchAddress(("Tag" + hName + "PluspxKalmanFit").c_str(), &m_hPlus_P_KalmanFit[0]);
  Tree->SetBranchAddress(("Tag" + hName + "PluspyKalmanFit").c_str(), &m_hPlus_P_KalmanFit[1]);
  Tree->SetBranchAddress(("Tag" + hName + "PluspzKalmanFit").c_str(), &m_hPlus_P_KalmanFit[2]);
  Tree->SetBranchAddress(("Tag" + hName + "PlusenergyKalmanFit").c_str(), &m_hPlus_P_KalmanFit[3]);
  Tree->SetBranchAddress(("Tag" + hName + "Minuspx").c_str(), &m_hMinus_P[0]);
  Tree->SetBranchAddress(("Tag" + hName + "Minuspy").c_str(), &m_hMinus_P[1]);
  Tree->SetBranchAddress(("Tag" + hName + "Minuspz").c_str(), &m_hMinus_P[2]);
  Tree->SetBranchAddress(("Tag" + hName + "Minusenergy").c_str(), &m_hMinus_P[3]);
  Tree->SetBranchAddress(("Tag" + hName + "MinuspxKalmanFit").c_str(), &m_hMinus_P_KalmanFit[0]);
  Tree->SetBranchAddress(("Tag" + hName + "MinuspyKalmanFit").c_str(), &m_hMinus_P_KalmanFit[1]);
  Tree->SetBranchAddress(("Tag" + hName + "MinuspzKalmanFit").c_str(), &m_hMinus_P_KalmanFit[2]);
  Tree->SetBranchAddress(("Tag" + hName + "MinusenergyKalmanFit").c_str(), &m_hMinus_P_KalmanFit[3]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "px").c_str(), &m_K0_P[0]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "py").c_str(), &m_K0_P[1]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "pz").c_str(), &m_K0_P[2]); // bug to be fixed in BOSS
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "energy").c_str(), &m_K0_P[3]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "pxKalmanFit").c_str(), &m_K0_P_KalmanFit[0]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "pyKalmanFit").c_str(), &m_K0_P_KalmanFit[1]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "pzKalmanFit").c_str(), &m_K0_P_KalmanFit[2]);
  Tree->SetBranchAddress(("Tag" + m_K0Mode + "energyKalmanFit").c_str(), &m_K0_P_KalmanFit[3]);
}

std::pair<int, int> KKpipi_vs_K0hh_PhaseSpace::Bin() const {
  int K0hhBin = GetK0hhBin();
  return K0hhBin > 0 ? std::make_pair(KKpipiBin(), K0hhBin) : std::make_pair(-KKpipiBin(), -K0hhBin);
}

std::pair<int, int> KKpipi_vs_K0hh_PhaseSpace::TrueBin() {
  int K0hhBin = GetTrueK0hhBin();
  return K0hhBin > 0 ? std::make_pair(TrueKKpipiBin(), K0hhBin) : std::make_pair(-TrueKKpipiBin(), -K0hhBin);
}

int KKpipi_vs_K0hh_PhaseSpace::GetK0hhBin() const {
  double M2Plus = m_KalmanFitSuccess == 1 ? (m_K0_P_KalmanFit + m_hPlus_P_KalmanFit).M2() : (m_K0_P + m_hPlus_P).M2();
  double M2Minus = m_KalmanFitSuccess == 1 ? (m_K0_P_KalmanFit + m_hMinus_P_KalmanFit).M2() : (m_K0_P + m_hMinus_P).M2();
  int Bin = DalitzUtilities::LookUpBinNumber(M2Plus, M2Minus, m_BinningScheme);
  if(Bin != 0) {
    return Bin;    
  } else {
    return DalitzUtilities::GetMappedK0hhBin(M2Plus, M2Minus, m_BinningScheme);

  }
}

int KKpipi_vs_K0hh_PhaseSpace::GetTrueK0hhBin() {
  // Get index of signal and tag D in truth information
  FindDIndex();
  // Get the end index of the tag D
  std::vector<int>::size_type TagEnd_index;
  if(m_TrueKinematics.SignalD_index < m_TrueKinematics.TagD_index) {
    TagEnd_index = m_TrueKinematics.NumberParticles;
  } else {
    TagEnd_index = m_TrueKinematics.SignalD_index;
  }
  // Indices of K0, hPlus, hMinus in truth information
  std::vector<int>::size_type K0_index, hPlus_index, hMinus_index;
  // Particle ID of h+ and h-
  int h_ID = m_hhMode == "pipi" ? 211 : 321;
  // Particle ID of K0
  int K0_ID;
  if(m_KStoKLBackground || m_K0Mode == "KS") {
    K0_ID = 310;
  } else {
    K0_ID = 130;
  }
  // begin() iterator of true IDs
  auto ID_begin = m_TrueKinematics.ParticleIDs.begin();
  // Get index of K0 in truth information
  auto K0_iter = std::find(ID_begin + m_TrueKinematics.TagD_index + 1, ID_begin + TagEnd_index, K0_ID);
  if(K0_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find K0 in truth information");
  } else {
    K0_index = K0_iter - ID_begin;
  }
  // Number of daughters to skip between 130/310 and the D daughters
  int SkipDaughters = 0;
  if(m_KStoKLBackground) {
    SkipDaughters = 7;
  } else if(K0_ID == 310) {
    SkipDaughters = 3;
  } else if(K0_ID == 130) {
    SkipDaughters = 1;
  }
  // Repeat for hPlus
  auto hPlus_iter = std::find(K0_iter + SkipDaughters, ID_begin + TagEnd_index, h_ID);
  if(hPlus_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find h+ in truth information");
  } else {
    hPlus_index = hPlus_iter - ID_begin;
  }
  // Repeat for hMinus
  auto hMinus_iter = std::find(K0_iter + SkipDaughters, ID_begin + TagEnd_index, -h_ID);
  if(hMinus_iter == ID_begin + TagEnd_index) {
    throw std::runtime_error("Cannot find h- in truth information");
  } else {
    hMinus_index = hMinus_iter - ID_begin;
  }
  // Calculate Dalitz coordinates
  TLorentzVector K0_P(m_TrueKinematics.TruePx[K0_index], m_TrueKinematics.TruePy[K0_index], m_TrueKinematics.TruePz[K0_index], m_TrueKinematics.TrueEnergy[K0_index]);
  TLorentzVector hPlus_P(m_TrueKinematics.TruePx[hPlus_index], m_TrueKinematics.TruePy[hPlus_index], m_TrueKinematics.TruePz[hPlus_index], m_TrueKinematics.TrueEnergy[hPlus_index]);
  TLorentzVector hMinus_P(m_TrueKinematics.TruePx[hMinus_index], m_TrueKinematics.TruePy[hMinus_index], m_TrueKinematics.TruePz[hMinus_index], m_TrueKinematics.TrueEnergy[hMinus_index]);
  double M2Plus = (K0_P + hPlus_P).M2();
  double M2Minus = (K0_P + hMinus_P).M2();
  // Get bin number
  int Bin = DalitzUtilities::LookUpBinNumber(M2Plus, M2Minus, m_BinningScheme);
  if(Bin != 0) {
    return Bin;
  } else {
    return DalitzUtilities::GetMappedK0hhBin(M2Plus, M2Minus, m_BinningScheme);
  }
}
