// Martin Duy Tat 25th November 2021

#include<string>
#include<stdexcept>
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooAddPdf.h"
#include"RooAbsReal.h"
#include"RooFormulaVar.h"
#include"RooSimultaneous.h"
#include"RooDataSet.h"
#include"TTree.h"
#include"BinnedFitModel.h"
#include"Settings.h"
#include"Utilities.h"
#include"Unique.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"

BinnedFitModel::BinnedFitModel(const Settings &settings,
			       RooRealVar *SignalMBC,
			       RooRealVar *TagMBC): m_SignalMBC(SignalMBC),
						    m_TagMBC(TagMBC),
						    m_Category(settings),
						    m_Settings(settings) {
  // First initialize the simultaneous fit with the correct categories for all bins
  auto CategoryVariable = m_Category.GetCategoryVariable();
  m_PDF = new RooSimultaneous(("Simultaneous_PDF_KKpipi_vs_" + settings.get("Mode")).c_str(), "", *CategoryVariable);
  // Set up all yield variables
  InitializeYields();
  // Set up the signal shape using signal MC
  InitializeSignalShape();
  // Then set up the combinatorial shape using an Argus function
  InitializeArgusShape();
  // Set up peaking backgrounds
  InitializePeakingBackgroundShapes();
  // Finally set up the simultaneous PDF for all bins
  InitializePDF();
}

BinnedFitModel::~BinnedFitModel() {
  if(m_PDF) {
    delete m_PDF;
  }
  for(auto &PeakingBackgroundShape : m_PeakingBackgroundShapes) {
    delete PeakingBackgroundShape.second;
  }
  m_PeakingBackgroundShapes.clear();
}

RooSimultaneous* BinnedFitModel::GetPDF() {
  return m_PDF;
}

void BinnedFitModel::InitializeYields() {
  RooArgList SignalYieldVars;
  int Bins = 0;
  for(const auto &CategoryString : m_Category.GetCategories()) {
    m_Yields.insert({CategoryString + "_CombinatorialYield", Unique::create<RooRealVar*>((CategoryString + "_CombinatorialYield").c_str(), "", 1.0, 0.0, 1000.0)});
    m_Yields.insert({CategoryString + "_SignalYield", Unique::create<RooRealVar*>((CategoryString + "_SignalYield").c_str(), "", 10.0, 0.0, 1000.0)});
    SignalYieldVars.add(*m_Yields.at(CategoryString + "_SignalYield"));
    Bins++;
  }
  std::string Formula("");
  for(int i = 0; i < Bins; i++) {
    if(i != 0) {
      Formula += "+";
    }
    Formula += "@";
    Formula += std::to_string(i);
  }
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    double BackgroundSignalRatio = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio");
    // For inclusive fit, quantum correlation is accounted for with a simple correction factor
    if(m_Settings.getB("Inclusive_fit") && m_Settings["MBC_Shape"].contains(Name + "_QuantumCorrelationFactor")) {
      BackgroundSignalRatio *= m_Settings["MBC_Shape"].getD(Name + "_QuantumCorrelationFactor");
    }
    for(const auto &CategoryString : m_Category.GetCategories()) {
      double BinYieldRatio = BackgroundSignalRatio*m_Settings["MBC_Shape"].getD(Name + "_" + CategoryString + "_frac_Yield");
      RooAbsReal *PeakingYield = Unique::create<RooFormulaVar*>((CategoryString + "_PeakingBackground" + std::to_string(i) + "Yield").c_str(), Form(("%f*(" + Formula + ")").c_str(), BinYieldRatio), SignalYieldVars);
      m_Yields.insert({CategoryString + "_PeakingBackground" + std::to_string(i) + "Yield", PeakingYield});
    }
  }
}

void BinnedFitModel::InitializeSignalShape() {
  m_Parameters.insert({"Mean1", Unique::create<RooRealVar*>((m_Settings.get("Mode") + "_DoubleTag_Mean1").c_str(), "", 0.0, -0.001, 0.001)});
  m_Parameters.insert({"Sigma1", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_Sigma1")});
  auto Resolution = Unique::create<RooGaussian*>("Gaussian1", "", *m_SignalMBC, *m_Parameters["Mean1"], *m_Parameters["Sigma1"]);
  TChain SignalMCChain(m_Settings.get("TreeName").c_str());
  std::string SignalMCFilename = m_Settings["Datasets_WithDeltaECuts"].get("SignalMC_DT");
  SignalMCFilename = Utilities::ReplaceString(SignalMCFilename, "TAG", m_Settings.get("Mode"));
  SignalMCChain.Add(SignalMCFilename.c_str());
  SignalMCChain.SetBranchStatus("*", 0);
  SignalMCChain.SetBranchStatus("SignalMBC", 1);
  SignalMCChain.SetBranchStatus("TagMBC", 1);
  TTree *ClonedMCChain = nullptr;
  if(m_Settings.getI("Events_in_MC") < 0 || m_Settings.getI("Events_in_MC") > SignalMCChain.GetEntries()) {
    ClonedMCChain = &SignalMCChain;
  } else {
    ClonedMCChain = SignalMCChain.CloneTree(m_Settings.getI("Events_in_MC"));
  }
  RooDataSet MCSignal("MCSignal", "", ClonedMCChain, RooArgList(*m_SignalMBC, *m_TagMBC)/*, "TagMBC > 1.86 && TagMBC < 1.87"*/);
  auto SignalShape = Unique::create<RooKeysPdf*>("SignalShape", "", *m_SignalMBC, MCSignal);
  m_SignalShapeConv = Unique::create<RooFFTConvPdf*>("SignalShapeConv", "", *m_SignalMBC, *SignalShape, *Resolution);
}

void BinnedFitModel::InitializeArgusShape() {
  m_Parameters.insert({"End", Unique::create<RooRealVar*>("End", "", 1.8865)});
  m_Parameters.insert({"c", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_c")});
  m_Argus = Unique::create<RooArgusBG*>("Argus", "", *m_SignalMBC, *m_Parameters["End"], *m_Parameters["c"]);
}

void BinnedFitModel::InitializePeakingBackgroundShapes() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    std::string PeakingShape = m_Settings["MBC_Shape"].get(Name + "_Shape");
    FitShape *PeakingPDF = nullptr;
    if(PeakingShape == "DoubleGaussian") {
      PeakingPDF = new DoubleGaussian_Shape(Name, m_Settings["MBC_Shape"], m_SignalMBC);
    } else if(PeakingShape == "DoubleCrystalBall") {
      PeakingPDF = new DoubleCrystalBall_Shape(Name, m_Settings["MBC_Shape"], m_SignalMBC);
    } else {
      throw std::invalid_argument("Unknown peaking background shape");
    }
    m_PeakingBackgroundShapes.insert({Name, PeakingPDF});
  }
}

RooAddPdf* BinnedFitModel::CreateBinPDF(const std::string &CategoryString) {
  RooArgList Shapes(*m_SignalShapeConv, *m_Argus);
  RooArgList Yields(*m_Yields.at(CategoryString + "_SignalYield"), *m_Yields.at(CategoryString + "_CombinatorialYield"));
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    Shapes.add(*m_PeakingBackgroundShapes.at(Name)->GetPDF());
    Yields.add(*m_Yields.at(CategoryString + "_PeakingBackground" + std::to_string(i) + "Yield"));
  }
  return Unique::create<RooAddPdf*>((CategoryString + "_PDF").c_str(), "", Shapes, Yields);
}

void BinnedFitModel::InitializePDF() {
  for(const auto &CategoryString : m_Category.GetCategories()) {
    RooAddPdf* BinPDF = CreateBinPDF(CategoryString);
    m_PDF->addPdf(*BinPDF, CategoryString.c_str());
  }
}

double BinnedFitModel::GetFractionInSignalRegion() const {
  using namespace RooFit;
  m_SignalMBC->setRange("SignalRange", 1.86, 1.87);
  return m_SignalShapeConv->createIntegral(*m_SignalMBC, NormSet(*m_SignalMBC), Range("SignalRange"))->getVal();
}
