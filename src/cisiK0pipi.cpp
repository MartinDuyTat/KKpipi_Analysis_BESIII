// Martin Duy Tat 14th April 2022

#include<stdexcept>
#include<vector>
#include<string>
#include<fstream>
#include"TMatrixT.h"
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooMultiVarGaussian.h"
#include"RooArgList.h"
#include"cisiK0pipi.h"
#include"Unique.h"
#include"Settings.h"

void cisiK0pipi::Initialise(const Settings &settings) {
  // First initialise Ki
  InitialiseKi(settings, "KSpipi");
  InitialiseKi(settings, "KLpipi");
  // Then initialise ci and si
  Initialisecisi(settings);
  // Finally set flag to true
  Initialised = true;
}

void cisiK0pipi::InitialiseKi(const Settings &settings, const std::string &TagMode) {
  for(int i = 1; i <= 8; i++) {
    std::string pName = TagMode + "_K_p" + std::to_string(i);
    std::string mName = TagMode + "_K_m" + std::to_string(i);
    double Ki = settings[TagMode + "_BinningScheme"]["Ki"].getD(pName);
    double Ki_err = settings[TagMode + "_BinningScheme"]["Ki"].getD(pName + "_err");
    double Kbari = settings[TagMode + "_BinningScheme"]["Ki"].getD(mName);
    double Kbari_err = settings[TagMode + "_BinningScheme"]["Ki"].getD(mName + "_err");
    auto Ki_Mean = Unique::create<RooRealVar*>((pName + "_Mean").c_str(), "", Ki);
    auto Ki_Sigma = Unique::create<RooRealVar*>((pName + "_Sigma").c_str(), "", Ki_err);
    auto Ki_Var = Unique::create<RooRealVar*>((pName + "_Var").c_str(), "", Ki, 0.0, 1.0);
    auto Ki_Gaussian = Unique::create<RooGaussian*>((pName + "_Gaussian").c_str(), "", *Ki_Var, *Ki_Mean, *Ki_Sigma);
    auto Kbari_Mean = Unique::create<RooRealVar*>((mName + "_Mean").c_str(), "", Kbari);
    auto Kbari_Sigma = Unique::create<RooRealVar*>((mName + "_Sigma").c_str(), "", Kbari_err);
    auto Kbari_Var = Unique::create<RooRealVar*>((mName + "_Var").c_str(), "", Kbari, 0.0, 1.0);
    auto Kbari_Gaussian = Unique::create<RooGaussian*>((mName + "_Gaussian").c_str(), "", *Kbari_Var, *Kbari_Mean, *Kbari_Sigma);
    Ki_Var->setConstant(true);
    Kbari_Var->setConstant(true);
    m_GaussianConstraintPDFs.push_back(Ki_Gaussian);
    m_GaussianConstraintPDFs.push_back(Kbari_Gaussian);
    if(TagMode == "KSpipi") {
      m_Ki_KSpipi.add(*Ki_Var);
      m_Ki_KSpipi.add(*Kbari_Var);
    } else if(TagMode == "KLpipi") {
      m_Ki_KLpipi.add(*Ki_Var);
      m_Ki_KLpipi.add(*Kbari_Var);
    } else {
      throw std::invalid_argument(TagMode + " tag mode is not KSpipi or KLpipi");
    }
  }
}

void cisiK0pipi::Initialisecisi(const Settings &settings) {
  RooArgList cisi_Mean;
  std::vector<double> cisi_Sigma;
  std::vector<std::string> TagModes{"KSpipi", "KLpipi"};
  std::vector<std::string> c_or_s{"c", "s"};
  for(const auto &Tag : TagModes) {
    for(const auto &cs : c_or_s) {
      for(int i = 1; i <= 8; i++) {
	std::string Name = Tag + "_" + cs + std::to_string(i);
	double cisi = settings[Tag + "_BinningScheme"]["cisi"].getD(Name);
	double cisi_err = settings[Tag + "_BinningScheme"]["cisi"].getD(Name + "_err");
	auto Mean = Unique::create<RooRealVar*>((Name + "_Mean").c_str(), "", cisi);
	cisi_Sigma.push_back(cisi_err);
	auto Var = Unique::create<RooRealVar*>((Name + "_Var").c_str(), "", cisi, -1.5, 1.5);
	Var->setConstant(true);
	cisi_Mean.add(*Mean);
	m_cisi.add(*Var);
      }
    }
  }
  auto CorrMatrix = ParseCorrelationMatrix(settings);
  TMatrixDSym CovMatrix(32);
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < 32; j++) {
      CovMatrix(i, j) = CorrMatrix(i, j)*cisi_Sigma[i]*cisi_Sigma[j];
    }
  }
  auto cisi_Gaussian = Unique::create<RooMultiVarGaussian*>("K0pipi_cisi_Gaussian", "", m_cisi, cisi_Mean, CovMatrix);
  m_GaussianConstraintPDFs.push_back(cisi_Gaussian);
}

TMatrixT<double> cisiK0pipi::ParseCorrelationMatrix(const Settings &settings) const {
  TMatrixT<double> CorrMatrix(32, 32);
  std::ifstream CorrFile(settings.get("KSpipi_KLpipi_cisi_Correlation"));
  for(int i = 0; i < 32; i++) {
    std::string Line;
    std::getline(CorrFile, Line);
    std::stringstream ss(Line);
    for(int j = 0; j < 32; j++) {
      if(j >= i) {
	double Value;
	ss >> Value;
	Value /= 100.0;
	CorrMatrix(i, j) = Value;
      } else {
	CorrMatrix(i, j) = CorrMatrix(j, i);
      }
    }
  }
  CorrFile.close();
  return CorrMatrix;
}
