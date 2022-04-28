// Martin Duy Tat 14th December 2021

#include<iostream>
#include<string>
#include<fstream>
#include<utility>
#include<stdexcept>
#include"TString.h"
#include"TMatrixTSym.h"
#include"TMatrixT.h"
#include"TMath.h"
#include"TFile.h"
#include"TTree.h"
#include"TRandom.h"
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
#include"CholeskySmearing.h"

FPlusFitter::FPlusFitter(const Settings &settings): m_Settings(settings),
						    m_FPlus_Model(m_Settings["FPlus_TagModes"].getD("KKpipi")),
						    m_KKpipi_BF_PDG(m_Settings["BranchingFractions"].getD("KKpipi")),
						    m_FPlus("FPlus", "", 0.5, 0.0, 3.0),
						    m_KKpipi_BF_CP("KKpipi_BF_CP", "", m_KKpipi_BF_PDG, 0.000, 0.005),
						    m_KKpipi_BF_KSpipi("KKpipi_BF_KSpipi", "", m_KKpipi_BF_PDG, 0.000, 0.005),
						    m_KKpipi_BF_KLpipi("KKpipi_BF_KLpipi", "", m_KKpipi_BF_PDG, 0.000, 0.005),
                                                    m_RunMinos(m_Settings.getB("RunMinos")) {
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
  auto Result = Model->fitTo(Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs), RooFit::Minos(m_RunMinos));
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
  auto Result = Model->fitTo(*Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs), RooFit::Minos(m_RunMinos));
  Result->Print("V");
  SaveFitResults(Result);
}

void FPlusFitter::DoManyToysOrFits(RooMultiVarGaussian *Model, const std::string RunMode) {
  int Seed = m_Settings.getI("Seed");
  if(RunMode == "ManyToys") {
    std::cout << "Run mode: Many toys\n";
    RooRandom::randomGenerator()->SetSeed(Seed);
  } else if(RunMode == "ManyFits") {
    std::cout << "Run mode: Many fits\n";
    gRandom->SetSeed(Seed);
  } else {
    return;
  }
  std::string Filename = m_Settings.get(RunMode + "OutputFilename");
  TFile OutputFile(Filename.c_str(), "RECREATE");
  TTree Tree("FPlusTree", "");
  int Status, CovQual;
  double FPlus, FPlus_err, FPlus_pull, Norm_CP, Norm_KSpipi, Norm_KLpipi, Norm_CP_pull, Norm_KSpipi_pull, Norm_KLpipi_pull;
  Tree.Branch("Status", &Status);
  Tree.Branch("CovQual", &CovQual);
  Tree.Branch("FPlus", &FPlus);
  Tree.Branch("FPlus_err", &FPlus_err);
  Tree.Branch("FPlus_pull", &FPlus_pull);
  Tree.Branch("KKpipi_BF_CP", &Norm_CP);
  Tree.Branch("KKpipi_BF_KSpipi", &Norm_KSpipi);
  Tree.Branch("KKpipi_BF_KLpipi", &Norm_KLpipi);
  Tree.Branch("KKpipi_BF_CP_pull", &Norm_CP_pull);
  Tree.Branch("KKpipi_BF_KSpipi_pull", &Norm_KSpipi_pull);
  Tree.Branch("KKpipi_BF_KLpipi_pull", &Norm_KLpipi_pull);
  int nToys = m_Settings.getI("NumberRuns");
  for(int i = 0; i < nToys; i++) {
    std::cout << "Run number " << i << "\n";
    ResetParameters();
    // Generate or smear dataset
    RooFitResult *Result = nullptr;
    if(RunMode == "ManyToys") {
      RooDataSet *Data = Model->generate(m_NormalizedYields, m_Settings.getI("StatsMultiplier"));
      Result = Model->fitTo(*Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs), RooFit::Minos(m_RunMinos));
      Data->Print("V");
    } else {
      ResetMeasurements();
      RooDataSet Data("Data", "", m_NormalizedYields);
      Data.add(m_NormalizedYields);
      Result = Model->fitTo(Data, RooFit::Save(), RooFit::ExternalConstraints(m_GaussianConstraintPDFs), RooFit::Minos(m_RunMinos));
      Data.Print("V");
    }
    Result->Print("V");
    Status = Result->status();
    CovQual = Result->covQual();
    FPlus = m_FPlus.getVal();
    FPlus_err = m_FPlus.getError();
    FPlus_pull = (FPlus - m_FPlus_Model)/FPlus_err;
    Norm_CP = m_KKpipi_BF_CP.getVal();
    Norm_KSpipi = m_KKpipi_BF_KSpipi.getVal();
    Norm_KLpipi = m_KKpipi_BF_KLpipi.getVal();
    Norm_CP_pull = (Norm_CP - m_KKpipi_BF_PDG)/m_KKpipi_BF_CP.getError();
    Norm_KSpipi_pull = (Norm_KSpipi - m_KKpipi_BF_PDG)/m_KKpipi_BF_KSpipi.getError();
    Norm_KLpipi_pull = (Norm_KLpipi - m_KKpipi_BF_PDG)/m_KKpipi_BF_KLpipi.getError();
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

void FPlusFitter::AddMeasurement_CP(const std::string &TagMode, bool Smearing) {
  // Get raw single tag yields and efficiency
  auto [ST_Yield, ST_Yield_err] = GetTagYield(TagMode, "ST", Smearing);
  double ST_Eff = GetEfficiency(TagMode, "ST", Smearing);
  ST_Yield /= ST_Eff;
  ST_Yield_err /= ST_Eff;
  // Get raw double tag yields and efficiency
  auto [DT_Yield, DT_Yield_err] = GetTagYield(TagMode, "DT", Smearing);
  double DT_Eff;
  if(m_Settings.get("Systematics") == "EfficiencyFactorisation") {
    DT_Eff = GetEfficiency("KKpipi", "ST", Smearing)*GetEfficiency(TagMode, "ST", Smearing);
  } else {
    DT_Eff = GetEfficiency(TagMode, "DT", Smearing);
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

void FPlusFitter::AddMeasurement_KShh(const std::string &TagMode, bool Smearing) {
  // Get raw single tag yields and efficiency
  auto [ST_Yield, ST_Yield_err] = GetTagYield(TagMode, "ST", Smearing);
  double ST_Eff = GetEfficiency(TagMode, "ST", Smearing);
  ST_Yield /= ST_Eff;
  ST_Yield_err /= ST_Eff;
  // Get number of bins
  int Bins = m_Settings[TagMode + "_BinningScheme"].getI("NumberBins");
  auto [DT_Yields_EffCorrected, DT_Yields_err] = GetBinnedTagYield(TagMode, Smearing);
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

void FPlusFitter::ResetMeasurements() {
  m_NormalizedYields.removeAll();
  m_Uncertainties.clear();
  for(const auto &Tag : m_TagModes) {
    if(Tag == "KSpipi" || Tag == "KSKK" || Tag == "KLpipi" || Tag == "KLKK" || Tag == "KSpipiPartReco") {
      AddMeasurement_KShh(Tag, true);
    } else {
      AddMeasurement_CP(Tag, true);
    }
  }
}

double FPlusFitter::GetEfficiency(std::string TagMode, const std::string &TagType, bool Smearing) const {
  std::string SingleDouble = TagType == "ST" ? "Single" : "Double";
  if(TagMode == "KSomega") {
    TagMode = "KSpipipi0";
  }
  double Eff = m_Settings[TagType + "_Efficiency"].getD(TagMode + "_" + SingleDouble + "TagEfficiency");
  if(Smearing && m_Settings.get("Systematics") == "Efficiency") {
    // Single tag efficiencies for KL modes cancel in the equation
    if(!(TagType == "ST" && TagMode.substr(0, 2) == "KL")) {
      double Eff_err = m_Settings[TagType + "_Efficiency"].getD(TagMode + "_" + SingleDouble + "TagEfficiency_err");
      Eff += gRandom->Gaus(0.0, Eff_err);
    }
  }
  return Eff;
}

TMatrixT<double>* FPlusFitter::GetEfficiencyMatrix(const std::string &TagMode, bool Smearing) const {
  TFile EffMatrixFile(m_Settings.get(TagMode + "_EfficiencyMatrix").c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  EffMatrixFile.GetObject("EffMatrix", EffMatrix);
  if(Smearing && m_Settings.get("Systematics") == "Efficiency") {
    TMatrixT<double> *EffMatrix_err = nullptr;
    EffMatrixFile.GetObject("EffMatrix_err", EffMatrix_err);
    for(int i = 0; i < EffMatrix->GetNrows(); i++) {
      for(int j = 0; j < EffMatrix->GetNcols(); j++) {
	(*EffMatrix)(i, j) += gRandom->Gaus(0.0, (*EffMatrix_err)(i, j));
      }
    }
  }
  EffMatrixFile.Close();
  EffMatrix->Invert();
  return EffMatrix;
}

std::pair<double, double> FPlusFitter::GetTagYield(const std::string &TagMode, const std::string &TagType, bool Smearing) const {
  std::string YieldName, SettingsName;
  if(TagType == "ST") {
    YieldName = TagMode + "_SingleTag_Yield";
    SettingsName = TagMode + "_ST_Yield";
  } else if(TagType == "DT") {
    YieldName = "DoubleTag_CP_KKpipi_vs_" + TagMode + "_SignalBin0_SignalYield";
    SettingsName = TagMode + "_DT_Yield";
  } else {
    throw std::invalid_argument(TagType + " is not a recognized tag type");
  }
  double Yield = m_Settings[SettingsName].getD(YieldName);
  double Yield_err;
  if(TagType == "ST" && TagMode.substr(0, 2) == "KL") {
    Yield_err = 0.0;
  } else { 
    Yield_err = m_Settings[SettingsName].getD(YieldName + "_err");
  }
  if(Smearing && m_Settings.get("Systematics") == "PeakingBackgrounds" && TagMode.substr(0, 2) != "KL") {
    double YieldSystError = m_Settings[SettingsName].getD(YieldName + "_PeakingBackgrounds_syst_err");
    Yield += gRandom->Gaus(0.0, YieldSystError);
  } else if(Smearing && m_Settings.get("Systematics") == "KL_ST_Yield" && TagMode.substr(0, 2) == "KL") {
    double YieldSystError = m_Settings[SettingsName].getD(YieldName + "_err");
    Yield += gRandom->Gaus(0.0, YieldSystError);
  }
  return std::make_pair(Yield, Yield_err);
}

std::pair<TMatrixT<double>, TMatrixT<double>> FPlusFitter::GetBinnedTagYield(const std::string &TagMode, bool Smearing) {
  int Bins = m_Settings[TagMode + "_BinningScheme"].getI("NumberBins");
  // Get double tag yield in each bin and correct for bin migration
  TMatrixT<double> DT_Yields(Bins, 1);
  TMatrixT<double> DT_Yields_err(Bins, 1);
  for(int i = 0; i < Bins; i++) {
    std::string DT_Name("DoubleTag_SCMB_KKpipi_vs_" + TagMode + "_SignalBin0_TagBin" + std::to_string(i + 1) + "_SignalYield");
    DT_Yields(i, 0) = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name);
    DT_Yields_err(i, 0) = m_Settings[TagMode + "_DT_Yield"].getD(DT_Name + "_err");
  }
  if(Smearing && m_Settings.get("Systematics") == "PeakingBackgrounds") {
    SmearBinnedTagYield(TagMode, DT_Yields);
  }
  auto EffMatrix = GetEfficiencyMatrix(TagMode, Smearing);
  TMatrixT<double> DT_Yields_EffCorrected = *EffMatrix*DT_Yields;
  for(int i = 0; i < Bins; i++) {
    DT_Yields_err(i, 0) *= (*EffMatrix)(i, i);
  }
  return std::make_pair(DT_Yields_EffCorrected, DT_Yields_err);
}

void FPlusFitter::SmearBinnedTagYield(const std::string &TagMode, TMatrixT<double> &DT_Yields) {
  if(m_CholeskySmearings.find(TagMode) == m_CholeskySmearings.end()) {
    std::string CovMatrixFilename = m_Settings.get(TagMode + "_PeakingBackground_CovMatrix");
    TFile CovMatrixFile(CovMatrixFilename.c_str(), "READ");
    TMatrixT<double> *CovMatrix = nullptr;
    CovMatrixFile.GetObject("CovMatrix", CovMatrix);
    m_CholeskySmearings.insert({TagMode, CholeskySmearing(*CovMatrix)});
  }
  m_CholeskySmearings.at(TagMode).Smear();
  DT_Yields += m_CholeskySmearings.at(TagMode).GetSmearings();
}
