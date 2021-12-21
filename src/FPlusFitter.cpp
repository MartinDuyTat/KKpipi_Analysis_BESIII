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
#include"RooFitResult.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"Settings.h"
#include"Unique.h"
#include"FPlusFitter.h"

FPlusFitter::FPlusFitter(const Settings &settings): m_FPlus("FPlus", "", 0.0, 1.0),
						    m_KKpipi_BF("KKpipi_BF", "", settings["BranchingFractions"].getD("KKpipi"), 0.000, 0.004),
						    m_Settings(settings) {
  if(!m_Settings.getB("Float_KKpipi_BF")) {
    m_KKpipi_BF.setConstant();
  }
}

void FPlusFitter::AddTag(const std::string &TagMode) {
  if(TagMode == "KSpipi" || TagMode == "KSKK") {
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
  auto Result = Model.fitTo(Data, RooFit::Save());
  Result->Print("V");
  SaveFitResults(Result);
}

void FPlusFitter::AddMeasurement_CP(const std::string &TagMode) {
  // Get raw single tag yields and efficiency
  double ST_Yield = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield");
  double ST_Yield_err = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield_err");
  // Get raw double tag yields and efficiency
  std::string DT_Name("DoubleTag_CP_KKpipi_vs_" + TagMode + "_SignalBin0_SignalYield");
  double DT_Yield = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name);
  double DT_Yield_err = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name + "_err");
  // Create variable with a normalized yield and add to dataset
  auto NormalizedYield = Unique::create<RooRealVar*>((TagMode + "_Normalized_Yield").c_str(), "", DT_Yield/ST_Yield);
  m_NormalizedYields.add(*NormalizedYield);
  // Save statistical uncertainty
  m_Uncertainties.push_back(TMath::Sqrt(TMath::Power(DT_Yield_err/DT_Yield, 2)
                                      + TMath::Power(ST_Yield_err/ST_Yield, 2))*DT_Yield/ST_Yield);
  // Print data point
  std::cout << "Adding " << TagMode << " tag mode\n";
  std::cout << "DT/ST ratio: (" << 1000.0*DT_Yield/ST_Yield << " \u00b1 " << 1000.0*m_Uncertainties.back() << ")e-3\n";
}

void FPlusFitter::AddPrediction_CP(const std::string &TagMode) {
  double FPlus_Tag = m_Settings["FPlus_TagModes"].getD(TagMode);
  TString Formula(Form("@1*(1 - (2*%f - 1)*(2*@0 - 1))", FPlus_Tag));
  auto PredictedYield = Unique::create<RooFormulaVar*>((TagMode + "_Normalized_Yield_Prediction").c_str(), Formula, RooArgList(m_FPlus, m_KKpipi_BF));
  m_PredictedYields.add(*PredictedYield);
}

void FPlusFitter::AddMeasurement_KShh(const std::string &TagMode) {
  // Get raw single tag yields and efficiency
  double ST_Yield = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield");
  double ST_Yield_err = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield_err");
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
  EffMatrixFile.Close();
  std::cout << "Adding " << TagMode << " tag mode\n";
  for(int i = 0; i < Bins; i++) {
    auto NormalizedYield = Unique::create<RooRealVar*>((TagMode + "_Normalized_Yield_Bin" + std::to_string(i + 1)).c_str(), "", DT_Yields_EffCorrected(i, 0)/ST_Yield);
    m_NormalizedYields.add(*NormalizedYield);
    m_Uncertainties.push_back(TMath::Sqrt(TMath::Power(DT_Yields_err(i, 0)/DT_Yields_EffCorrected(i, 0), 2)
                                        + TMath::Power(ST_Yield_err/ST_Yield, 2))*DT_Yields_EffCorrected(i, 0)/ST_Yield);
    std::cout << "DT/ST ratio in bin " << i + 1 << ": (" << 1000.0*DT_Yields_EffCorrected(i, 0)/ST_Yield << " \u00b1 " << 1000.0*m_Uncertainties.back() << ")e-3\n";
  }
}

void FPlusFitter::AddPrediction_KShh(const std::string &TagMode) {
  int Bins = m_Settings[TagMode + "_BinningScheme"].getI("NumberBins");
  for(int i = 1; i <= Bins; i++) {
    double ci = m_Settings[TagMode + "_BinningScheme"]["cisi"].getD(TagMode + "_c" + std::to_string(i));
    double Ki = m_Settings[TagMode + "_BinningScheme"]["Ki"].getD(TagMode + "_K_p" + std::to_string(i));
    double Kbari = m_Settings[TagMode + "_BinningScheme"]["Ki"].getD(TagMode + "_K_m" + std::to_string(i));
    TString Formula(Form("0.5*@1*(%f + %f - 2*%f*sqrt(%f*%f)*(2*@0 - 1))", Ki, Kbari, ci, Ki, Kbari));
    auto PredictedYield = Unique::create<RooFormulaVar*>((TagMode + "_Normalized_Yield_Prediction_Bin" + std::to_string(i)).c_str(), Formula, RooArgList(m_FPlus, m_KKpipi_BF));
    m_PredictedYields.add(*PredictedYield);
  }
}

void FPlusFitter::SaveFitResults(RooFitResult *Result) const {
  std::ofstream Outfile(m_Settings.get("ResultsFile"));
  Outfile << "status " << Result->status() << "\n";
  Outfile << "covQual " << Result->covQual() << "\n\n";
  Outfile << "FPlus         " << m_FPlus.getVal() << "\n";
  Outfile << "FPlus_err     " << m_FPlus.getError() << "\n";
  Outfile << "BF_KKpipi     " << m_KKpipi_BF.getVal() << "\n";
  Outfile << "BF_KKpipi_err " << m_KKpipi_BF.getError() << "\n";
  Outfile.close();
}
