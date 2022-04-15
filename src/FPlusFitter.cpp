// Martin Duy Tat 14th December 2021

#include<iostream>
#include<string>
#include<fstream>
#include"TString.h"
#include"TMatrixTSym.h"
#include"TMatrixT.h"
#include"TMath.h"
#include"TFile.h"
#include"TTree.h"
#include"RooRealVar.h"
#include"RooFormulaVar.h"
#include"RooArgList.h"
#include"RooMultiVarGaussian.h"
#include"RooGaussian.h"
#include"RooFitResult.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooRandom.h"
#include"Settings.h"
#include"Unique.h"
#include"FPlusFitter.h"

FPlusFitter::FPlusFitter(const Settings &settings): m_Settings(settings),
						    m_FPlus_Model(m_Settings["FPlus_TagModes"].getD("KKpipi")),
						    m_KKpipi_BF_PDG(m_Settings["BranchingFractions"].getD("KKpipi")),
						    m_FPlus("FPlus", "", 0.5, 0.0, 3.0),
						    m_KKpipi_BF_CP("KKpipi_BF_CP", "", m_KKpipi_BF_PDG, 0.000, 0.005),
						    m_KKpipi_BF_KSpipi("KKpipi_BF_KSpipi", "", m_KKpipi_BF_PDG, 0.000, 0.005),
						    m_KKpipi_BF_KLpipi("KKpipi_BF_KLpipi", "", m_KKpipi_BF_PDG, 0.000, 0.005) {
  m_KKpipi_BF_CP.setConstant(true);
  m_KKpipi_BF_KSpipi.setConstant(true);
  m_KKpipi_BF_KLpipi.setConstant(true);
}

void FPlusFitter::AddTag(const std::string &TagMode) {
  m_TagModes.push_back(TagMode);
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
  std::string RunMode = m_Settings.get("RunMode");
  if(RunMode == "SingleFit") {
    DoSingleFit(&Model);
  } else if(RunMode == "SingleToy") {
    DoSingleToy(&Model);
  } else if(RunMode == "ManyToys") {
    DoManyToysOrFits(&Model, RunMode);
  } else if(RunMode == "ManyFits") {
    DoManyToysOrFits(&Model, RunMode);
  }
}

void FPlusFitter::DoSingleFit(RooMultiVarGaussian *Model) {
  std::cout << "Run mode: Single fit\n";
  // Set up dataset
  RooDataSet Data("Data", "", m_NormalizedYields);
  Data.add(m_NormalizedYields);
  Data.Print("V");
  auto Result = Model->fitTo(Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs));
  Result->Print("V");
  SaveFitResults(Result);
}

void FPlusFitter::DoSingleToy(RooMultiVarGaussian *Model) {
  std::cout << "Run mode: Single toy\n";
  // Set random seed
  int Seed = m_Settings.getI("ToySeed");
  RooRandom::randomGenerator()->SetSeed(Seed);
  // Generate dataset
  ResetParameters();
  RooDataSet *Data = Model->generate(m_NormalizedYields, m_Settings.getI("StatsMultiplier"));
  Data->Print("V");
  auto Result = Model->fitTo(*Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs));
  Result->Print("V");
  SaveFitResults(Result);
}

void FPlusFitter::DoManyToysOrFits(RooMultiVarGaussian *Model, const std::string RunMode) {
  if(RunMode == "ManyToys") {
    std::cout << "Run mode: Many toys\n";
  } else if(RunMode == "ManyFits") {
    std::cout << "Run mode: Many fits\n";
  } else {
    return;
  }
  std::string Filename = m_Settings.get(RunMode + "OutputFilename");
  TFile OutputFile(Filename.c_str(), "RECREATE");
  TTree Tree("FPlusTree", "");
  int Status, CovQual;
  double FPlus, FPlus_err, FPlus_pull;
  Tree.Branch("Status", &Status);
  Tree.Branch("CovQual", &CovQual);
  Tree.Branch("FPlus", &FPlus);
  Tree.Branch("FPlus_err", &FPlus_err);
  Tree.Branch("FPlus_pull", &FPlus_pull);
  int Seed = m_Settings.getI("Seed");
  int nToys = m_Settings.getI("NumberRuns");
  for(int i = 0; i < nToys; i++) {
    std::cout << "Run number " << i << "\n";
    ResetParameters();
    // Generate or smear dataset
    RooDataSet Data;
    RooFitResult *Result = nullptr;
    if(RunMode == "ManyToys") {
      // Set random seed
      RooRandom::randomGenerator()->SetSeed(Seed + i);
      Data = *Model->generate(m_NormalizedYields, m_Settings.getI("StatsMultiplier"));
    } else {
      ResetMeasurements(Seed + i);
      Data = RooDataSet("Data", "", m_NormalizedYields);
      Data.add(m_NormalizedYields);
    }
    Data.Print("V");
    Result = Model->fitTo(Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs));
    Result->Print("V");
    Status = Result->status();
    CovQual = Result->covQual();
    FPlus = m_FPlus.getVal();
    FPlus_err = m_FPlus.getError();
    FPlus_pull = (FPlus - m_FPlus_Model)/FPlus_err;
    Tree.Fill();
  }
  OutputFile.cd();
  Tree.Write();
  OutputFile.Close();
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
  std::string YieldName = TagMode + "_Normalized_Yield";
  m_YieldVars[YieldName] = RooRealVar(YieldName.c_str(), "", DT_Yield/ST_Yield);
  if(m_YieldVars[YieldName].getVal() < 0.0) {
    m_YieldVars[YieldName].setVal(0.0);
  }
  m_NormalizedYields.add(m_YieldVars[YieldName]);
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
  std::string YieldName = TagMode + "_Normalized_Yield_Prediction";
  RooArgList ParameterList(m_FPlus, m_KKpipi_BF_CP, *FPlus_Tag);
  auto PredictedYield = Unique::create<RooFormulaVar*>(YieldName.c_str(), Formula, ParameterList);
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
    std::string YieldName = TagMode + "_Normalized_Yield_Bin" + std::to_string(i + 1);
    m_YieldVars[YieldName] = RooRealVar(YieldName.c_str(), "", DT_Yields_EffCorrected(i, 0)/ST_Yield);
    if(m_YieldVars[YieldName].getVal() < 0.0) {
      m_YieldVars[YieldName].setVal(0.0);
    }
    m_NormalizedYields.add(m_YieldVars[YieldName]);
    double Uncertainty2 = TMath::Power(DT_Yields_err(i, 0)/DT_Yields_EffCorrected(i, 0), 2);
    Uncertainty2 += TMath::Power(ST_Yield_err/ST_Yield, 2);
    m_Uncertainties.push_back(TMath::Sqrt(Uncertainty2)*DT_Yields_EffCorrected(i, 0)/ST_Yield);
    std::cout << "DT/ST ratio in bin " << i + 1 << ": (";
    std::cout << 1000.0*DT_Yields_EffCorrected(i, 0)/ST_Yield << " \u00b1 " << 1000.0*m_Uncertainties.back() << ")e-3\n";
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
  auto FPlus_Tag = GetFPlusTag(TagMode == "KLpipi" ? "KLpipi" : "KSpipi");
  for(int i = 1; i <= Bins; i++) {
    std::string FormulaString;
    if(TagMode.substr(0, 2) == "KS") {
      FormulaString = "@1*(@2 + @3 - 2*@4*sqrt(@2*@3)*(2*@0 - 1))/(1 - (2*@5 - 1)*%f)";
    } else {
      FormulaString = "@1*(@2 + @3 + 2*@4*sqrt(@2*@3)*(2*@0 - 1))/(1 - (2*@5 - 1)*%f)";
    }
    double y_CP = m_Settings.getD("y_CP");
    TString Formula(Form(FormulaString.c_str(), y_CP));
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
    K0hhParameters.add(*FPlus_Tag);
    if(m_Settings.getB("GaussianConstrainExternalParameters")) {
      static_cast<RooRealVar*>(K0hhParameters.at(2))->setConstant(false);
      static_cast<RooRealVar*>(K0hhParameters.at(3))->setConstant(false);
      static_cast<RooRealVar*>(K0hhParameters.at(4))->setConstant(false);
      static_cast<RooRealVar*>(K0hhParameters.at(5))->setConstant(false);
    } else {
      static_cast<RooRealVar*>(K0hhParameters.at(5))->setConstant(true);
    }
    std::string YieldName = TagMode + "_Normalized_Yield_Prediction_Bin" + std::to_string(i);
    auto PredictedYield = Unique::create<RooFormulaVar*>(YieldName.c_str(), Formula, K0hhParameters);
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

void FPlusFitter::ResetParameters() {
  m_FPlus.setVal(m_FPlus_Model);
  m_KKpipi_BF_CP.setVal(m_KKpipi_BF_PDG);
  m_KKpipi_BF_KSpipi.setVal(m_KKpipi_BF_PDG);
  m_KKpipi_BF_KLpipi.setVal(m_KKpipi_BF_PDG);
}

void FPlusFitter::ResetMeasurements(int) {
  m_NormalizedYields.removeAll();
  m_Uncertainties.clear();
  for(const auto &Tag : m_TagModes) {
    if(Tag == "KSpipi" || Tag == "KSKK" || Tag == "KLpipi" || Tag == "KLKK" || Tag == "KSpipiPartReco") {
      AddMeasurement_KShh(Tag);
    } else {
      AddMeasurement_CP(Tag);
    }
  }
}
