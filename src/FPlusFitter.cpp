// Martin Duy Tat 14th December 2021

#include<iostream>
#include<string>
#include<fstream>
#include"TString.h"
#include"TMatrixTSym.h"
#include"TMatrixT.h"
#include"TMath.h"
#include"TFile.h"
#include"RooRealVar.h"
#include"RooFormulaVar.h"
#include"RooArgList.h"
#include"RooMultiVarGaussian.h"
#include"RooGaussian.h"
#include"RooFitResult.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"Settings.h"
#include"Unique.h"
#include"FPlusFitter.h"

FPlusFitter::FPlusFitter(const Settings &settings): m_FPlus("FPlus", "", 0.0, 1.0),
						    m_KKpipi_BF_CP("KKpipi_BF_CP", "", settings["BranchingFractions"].getD("KKpipi"), 0.000, 0.005),
						    m_KKpipi_BF_KSpipi("KKpipi_BF_KSpipi", "", settings["BranchingFractions"].getD("KKpipi"), 0.000, 0.005),
						    m_KKpipi_BF_KLpipi("KKpipi_BF_KLpipi", "", settings["BranchingFractions"].getD("KKpipi"), 0.000, 0.005),
						    m_Settings(settings) {
  m_KKpipi_BF_CP.setConstant(true);
  m_KKpipi_BF_KSpipi.setConstant(true);
  m_KKpipi_BF_KLpipi.setConstant(true);
}

void FPlusFitter::AddTag(const std::string &TagMode) {
  if(TagMode == "KSpipi" || TagMode == "KSKK" || TagMode == "KLpipi" || TagMode == "KLKK" || TagMode == "KSpipiPartReco") {
    AddMeasurement_KShh(TagMode);
    AddPrediction_KShh(TagMode);
  } else {
    AddMeasurement_CP(TagMode);
    AddPrediction_CP(TagMode);
  }
}

void FPlusFitter::InitializeAndFit() {
  // Set up trivial (diagonal, assuming all tags are independent) covariance matrix
  auto NumberTags = m_Uncertainties.size();
  TMatrixTSym<double> CovMatrix(NumberTags);
  for(std::size_t i = 0; i < NumberTags; i++) {
    CovMatrix(i, i) = m_Uncertainties[i]*m_Uncertainties[i];
  }
  // Set up multidimensional Gaussian containing predicted yields and normalized yields
  RooMultiVarGaussian Model("Model", "", m_NormalizedYields, m_PredictedYields, CovMatrix);
  // Set up dataset
  RooDataSet Data("Data", "", m_NormalizedYields);
  Data.add(m_NormalizedYields);
  Data.Print("V");
  // Perform fit
  auto Result = Model.fitTo(Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs));
  Result->Print("V");
  SaveFitResults(Result);
}

RooRealVar* FPlusFitter::GetFPlusTag(const std::string &TagMode) {
  double FPlus_Tag = m_Settings["FPlus_TagModes"].getD(TagMode);
  auto Mean = Unique::create<RooRealVar*>((TagMode + "_FPlus_Mean").c_str(), "", FPlus_Tag);
  RooRealVar *Sigma = nullptr;
  if(m_Settings["FPlus_TagModes"].contains(TagMode + "_err")) {
    double FPlus_Tag_Uncertainty = m_Settings["FPlus_TagModes"].getD(TagMode + "_err");
    Sigma = Unique::create<RooRealVar*>((TagMode + "_FPlus_Sigma").c_str(), "", FPlus_Tag_Uncertainty);
    auto FPlusVar = Unique::create<RooRealVar*>((TagMode + "_FPlus_Var").c_str(), "", FPlus_Tag, 0.0, 1.0);
    auto FPlusGaussian = Unique::create<RooGaussian*>((TagMode + "_Gaussian").c_str(), "", *FPlusVar, *Mean, *Sigma);
    m_GaussianConstraintPDFs.add(*FPlusGaussian);
    if(!m_Settings.getB("GaussianConstrainExternalParameters")) {
      FPlusVar->setConstant(true);
    } else {
      std::cout << "F+ of " << TagMode << " tag mode is Gaussian constrained to ";
      std::cout << FPlus_Tag << " \u00b1 " << FPlus_Tag_Uncertainty << "\n";
    }
    return FPlusVar;
  } else {
    return Mean;
  }
}

void FPlusFitter::AddMeasurement_CP(const std::string &TagMode) {
  // Get raw single tag yields and efficiency
  double ST_Yield = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield");
  double ST_Yield_err = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield_err");
  double ST_Eff;
  if(TagMode == "KSomega") {
    ST_Eff = m_Settings["ST_Efficiency"].getD("KSpipipi0_SingleTagEfficiency");
  } else {
    ST_Eff = m_Settings["ST_Efficiency"].getD(TagMode + "_SingleTagEfficiency");
  }
  ST_Yield /= ST_Eff;
  ST_Yield_err /= ST_Eff;
  // Get raw double tag yields and efficiency
  std::string DT_Name("DoubleTag_CP_KKpipi_vs_" + TagMode + "_SignalBin0_SignalYield");
  double DT_Yield = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name);
  double DT_Yield_err = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name + "_err");
  double DT_Eff;
  if(TagMode == "KSomega") {
    DT_Eff = m_Settings["DT_Efficiency"].getD("KSpipipi0_DoubleTagEfficiency");
  } else {
    DT_Eff = m_Settings["DT_Efficiency"].getD(TagMode + "_DoubleTagEfficiency");
  }
  DT_Yield /= DT_Eff;
  DT_Yield_err /= DT_Eff;
  // Create variable with a normalized yield and add to dataset
  auto NormalizedYield = Unique::create<RooRealVar*>((TagMode + "_Normalized_Yield").c_str(), "", DT_Yield/ST_Yield);
  if(NormalizedYield->getVal() < 0.0) {
    NormalizedYield->setVal(0.0);
  }
  m_NormalizedYields.add(*NormalizedYield);
  // Save statistical uncertainty
  m_Uncertainties.push_back(TMath::Sqrt(TMath::Power(DT_Yield_err/DT_Yield, 2)
                                      + TMath::Power(ST_Yield_err/ST_Yield, 2))*DT_Yield/ST_Yield);
  // Print data point
  std::cout << "Adding " << TagMode << " tag mode\n";
  std::cout << "DT/ST ratio: (" << 1000.0*DT_Yield/ST_Yield << " \u00b1 " << 1000.0*m_Uncertainties.back() << ")e-3\n";
}

void FPlusFitter::AddPrediction_CP(const std::string &TagMode) {
  if(m_Settings.getB("Float_KKpipi_BF")) {
    m_KKpipi_BF_CP.setConstant(false);
  }
  auto FPlus_Tag = GetFPlusTag(TagMode);
  double y_CP = m_Settings.getD("y_CP");
  TString Formula(Form("@1*(1 - (2*@2 - 1)*(2*@0 - 1))/(1 - (2*@2 - 1)*%f)", y_CP));
  auto PredictedYield = Unique::create<RooFormulaVar*>((TagMode + "_Normalized_Yield_Prediction").c_str(), Formula, RooArgList(m_FPlus, m_KKpipi_BF_CP, *FPlus_Tag));
  m_PredictedYields.add(*PredictedYield);
}

void FPlusFitter::AddMeasurement_KShh(const std::string &TagMode) {
  // Get raw single tag yields and efficiency
  double ST_Yield = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield");
  double ST_Yield_err = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield_err");
  double ST_Eff = m_Settings["ST_Efficiency"].getD(TagMode + "_SingleTagEfficiency");
  ST_Yield /= ST_Eff;
  ST_Yield_err /= ST_Eff;
  // Get number of bins
  int Bins = m_Settings[TagMode + "_BinningScheme"].getI("NumberBins");
  // Get double tag yield in each bin and correct for bin migration
  TMatrixT<double> DT_Yields(Bins, 1);
  TMatrixT<double> DT_Yields_err(Bins, 1);
  for(int i = 0; i < Bins; i++) {
    std::string DT_Name("DoubleTag_SCMB_KKpipi_vs_" + TagMode + "_SignalBin0_TagBin" + std::to_string(i + 1) + "_SignalYield");
    DT_Yields(i, 0) = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name);
    DT_Yields_err(i, 0) = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name + "_err");
  }
  TFile EffMatrixFile(m_Settings.get(TagMode + "_EfficiencyMatrix").c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  EffMatrixFile.GetObject("EffMatrix", EffMatrix);
  EffMatrix->Invert();
  TMatrixT<double> DT_Yields_EffCorrected = *EffMatrix*DT_Yields;
  for(int i = 0; i < Bins; i++) {
    DT_Yields_err(i, 0) *= (*EffMatrix)(i, i);
  }
  EffMatrixFile.Close();
  std::cout << "Adding " << TagMode << " tag mode\n";
  for(int i = 0; i < Bins; i++) {
    auto NormalizedYield = Unique::create<RooRealVar*>((TagMode + "_Normalized_Yield_Bin" + std::to_string(i + 1)).c_str(), "", DT_Yields_EffCorrected(i, 0)/ST_Yield);
    if(NormalizedYield->getVal() < 0.0) {
      NormalizedYield->setVal(0.0);
    }
    m_NormalizedYields.add(*NormalizedYield);
    m_Uncertainties.push_back(TMath::Sqrt(TMath::Power(DT_Yields_err(i, 0)/DT_Yields_EffCorrected(i, 0), 2)
                                        + TMath::Power(ST_Yield_err/ST_Yield, 2))*DT_Yields_EffCorrected(i, 0)/ST_Yield);
    std::cout << "DT/ST ratio in bin " << i + 1 << ": (" << 1000.0*DT_Yields_EffCorrected(i, 0)/ST_Yield << " \u00b1 " << 1000.0*m_Uncertainties.back() << ")e-3\n";
  }
}

void FPlusFitter::AddPrediction_KShh(const std::string &TagMode) {
  if(!m_cisi_K0pipi.Initialised) {
    m_cisi_K0pipi.Initialise(m_Settings);
  }
  for(auto PDF : m_cisi_K0pipi.m_GaussianConstraintPDFs) {
    m_GaussianConstraintPDFs.add(*PDF);
  }
  if(m_Settings.getB("Float_KKpipi_BF")) {
    if(TagMode == "KSpipi" || TagMode == "KSpipiPartReco") {
      m_KKpipi_BF_KSpipi.setConstant(false);
    } else if(TagMode == "KLpipi") {
      m_KKpipi_BF_KLpipi.setConstant(false);
    }
  }
  int Bins = m_Settings[TagMode + "_BinningScheme"].getI("NumberBins");
  std::string BinningTag = TagMode == "KSpipiPartReco" ? "KSpipi" : TagMode;
  for(int i = 1; i <= Bins; i++) {
    std::string Formula;
    if(TagMode.substr(0, 2) == "KS") {
      Formula = "@1*(@2 + @3 - 2*@4*sqrt(@2*@3)*(2*@0 - 1))";
    } else {
      Formula = "@1*(@2 + @3 + 2*@4*sqrt(@2*@3)*(2*@0 - 1))";
    }
    RooArgList K0hhParameters;
    K0hhParameters.add(m_FPlus);
    if(TagMode == "KLpipi") {
      K0hhParameters.add(m_KKpipi_BF_KLpipi);
      K0hhParameters.add(m_cisi_K0pipi.m_Ki_KLpipi[2*(i - 1)]); // Ki
      K0hhParameters.add(m_cisi_K0pipi.m_Ki_KLpipi[2*(i - 1) + 1]); // Kbari
      K0hhParameters.add(m_cisi_K0pipi.m_cisi[i - 1 + 16]); // ci
    } else {
      K0hhParameters.add(m_KKpipi_BF_KSpipi);
      K0hhParameters.add(m_cisi_K0pipi.m_Ki_KSpipi[2*(i - 1)]); // Ki
      K0hhParameters.add(m_cisi_K0pipi.m_Ki_KSpipi[2*(i - 1) + 1]); // Kbari
      K0hhParameters.add(m_cisi_K0pipi.m_cisi[i - 1]); // ci
    }
    if(m_Settings.getB("GaussianConstrainExternalParameters")) {
      static_cast<RooRealVar*>(K0hhParameters.at(2))->setConstant(false);
      static_cast<RooRealVar*>(K0hhParameters.at(3))->setConstant(false);
      static_cast<RooRealVar*>(K0hhParameters.at(4))->setConstant(false);
    }
    std::string YieldName = TagMode + "_Normalized_Yield_Prediction_Bin" + std::to_string(i);
    auto PredictedYield = Unique::create<RooFormulaVar*>(YieldName.c_str(), Formula.c_str(), K0hhParameters);
    m_PredictedYields.add(*PredictedYield);
  }
}

void FPlusFitter::SaveFitResults(RooFitResult *Result) const {
  std::ofstream Outfile(m_Settings.get("ResultsFile"));
  Outfile << "status               " << Result->status() << "\n";
  Outfile << "covQual              " << Result->covQual() << "\n\n";
  Outfile << "FPlus                " << m_FPlus.getVal() << "\n";
  Outfile << "FPlus_err            " << m_FPlus.getError() << "\n";
  if(!m_KKpipi_BF_CP.isConstant()) {
    Outfile << "BF_KKpipi_CP         " << m_KKpipi_BF_CP.getVal() << "\n";
    Outfile << "BF_KKpipi_CP_err     " << m_KKpipi_BF_CP.getError() << "\n";
    Outfile << "Correlation_CP       " << Result->correlation("FPlus", "KKpipi_BF_CP") << "\n";
  }
  if(!m_KKpipi_BF_KSpipi.isConstant()) {
    Outfile << "BF_KKpipi_KSpipi     " << m_KKpipi_BF_KSpipi.getVal() << "\n";
    Outfile << "BF_KKpipi_KSpipi_err " << m_KKpipi_BF_KSpipi.getError() << "\n";
    Outfile << "Correlation_KSpipi   " << Result->correlation("FPlus", "KKpipi_BF_KSpipi") << "\n";
  }
  if(!m_KKpipi_BF_KLpipi.isConstant()) {
    Outfile << "BF_KKpipi_KLpipi     " << m_KKpipi_BF_KLpipi.getVal() << "\n";
    Outfile << "BF_KKpipi_KLpipi_err " << m_KKpipi_BF_KLpipi.getError() << "\n";
    Outfile << "Correlation_KLpipi   " << Result->correlation("FPlus", "KKpipi_BF_KLpipi") << "\n";
  }
  Outfile << "MinLL                " << Result->minNll() << "\n";
  Outfile.close();
}
