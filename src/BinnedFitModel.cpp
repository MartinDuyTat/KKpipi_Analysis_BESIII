// Martin Duy Tat 25th November 2021

#include<string>
#include<stdexcept>
#include"TTree.h"
#include"TRandom.h"
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooAddPdf.h"
#include"RooAbsReal.h"
#include"RooFormulaVar.h"
#include"RooSimultaneous.h"
#include"RooDataSet.h"
#include"BinnedFitModel.h"
#include"Settings.h"
#include"Utilities.h"
#include"Unique.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"
#include"RooShapes/CrystalBall_Shape.h"
#include"RooShapes/Chebychev_Shape.h"

BinnedFitModel::BinnedFitModel(const Settings &settings,
			       RooRealVar *SignalMBC): m_SignalMBC(SignalMBC),
						       m_Category(settings),
						       m_Settings(settings) {
  // First initialize the simultaneous fit with the correct categories for all bins
  auto CategoryVariable = m_Category.GetCategoryVariable();
  m_PDF = new RooSimultaneous(("Simultaneous_PDF_KKpipi_vs_" + settings.get("Mode")).c_str(), "", *CategoryVariable);
  // Set up all yield variables
  InitializeYields();
  // Set up the signal shape using signal MC
  InitializeSignalShape();
  // Then set up the combinatorial shape
  InitializeCombinatorialShape();
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
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(const auto &CategoryString : m_Category.GetCategories()) {
    m_Yields.insert({CategoryString + "_CombinatorialYield", Unique::create<RooRealVar*>((CategoryString + "_CombinatorialYield").c_str(), "", 1.0, 0.0, 1000.0)});
    m_Yields.insert({CategoryString + "_SignalYield", Unique::create<RooRealVar*>((CategoryString + "_SignalYield").c_str(), "", 10.0, 0.0, 1000.0)});
    for(int i = 0; i < PeakingBackgrounds; i++) {
      std::string Name = Mode + "_PeakingBackground" + std::to_string(i) + "_" + CategoryString;
      if(!m_Settings["MBC_Shape"].contains(Name + "_Yield")) {
	auto BackgroundSignalRatio = Utilities::load_param(m_Settings["MBC_Shape"], Name + "_BackgroundToSignalRatio");
	std::cout << "Adding peaking background with background-to-signal ratio: " << BackgroundSignalRatio->getVal() << "\n";
	// Quantum correlation is accounted for with a correction factor
	RooRealVar *QCFactor = nullptr;
	if(m_Settings["MBC_Shape"].contains(Name + "_QuantumCorrelationFactor")) {
	  QCFactor = Utilities::load_param(m_Settings["MBC_Shape"], Name + "_QuantumCorrelationFactor");
	  std::cout << "Quantum correlation correction factor: " << QCFactor->getVal() << "\n";
	} else {
	  QCFactor = Unique::create<RooRealVar*>((Name + "QuantumCorrelationFactor").c_str(), "", 1.0);
	  std::cout << "Quantum correlation correction factor: None\n";
	}
	RooArgList PeakingYieldParameters(*m_Yields.at(CategoryString + "_SignalYield"), *BackgroundSignalRatio, *QCFactor);
	RooAbsReal *PeakingYield = Unique::create<RooFormulaVar*>((Name + "_Yield").c_str(), "@0*@1*@2", PeakingYieldParameters);
	m_Yields.insert({CategoryString + "_PeakingBackground" + std::to_string(i) + "Yield", PeakingYield});
      } else {
	double BackgroundYield = m_Settings["MBC_Shape"].getD(Name + "_Yield");
	RooAbsReal *PeakingYield = Unique::create<RooRealVar*>((Name + "_Yield").c_str(), "", BackgroundYield);
	std::cout << "Adding peaking background with yield: " << BackgroundYield << "\n";
	m_Yields.insert({CategoryString + "_PeakingBackground" + std::to_string(i) + "Yield", PeakingYield});
      }
    }
  }
}

void BinnedFitModel::InitializeSignalShape() {
  m_Parameters.insert({"Mean1", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_Mean1")});
  m_Parameters.insert({"Sigma1", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_Sigma1")});
  auto Resolution = Unique::create<RooGaussian*>("Gaussian1", "", *m_SignalMBC, *m_Parameters["Mean1"], *m_Parameters["Sigma1"]);
  TChain SignalMCChain(m_Settings.get("TreeName").c_str());
  std::string SignalMCFilename = m_Settings["Datasets_WithDeltaECuts"].get("SignalMC_DT");
  SignalMCFilename = Utilities::ReplaceString(SignalMCFilename, "TAG", m_Settings.get("Mode"));
  SignalMCChain.Add(SignalMCFilename.c_str());
  SignalMCChain.SetBranchStatus("*", 0);
  SignalMCChain.SetBranchStatus(m_Settings.get("FitVariable").c_str(), 1);
  TTree *ClonedMCChain = nullptr;
  if(m_Settings.getI("Events_in_MC") < 0 || m_Settings.getI("Events_in_MC") > SignalMCChain.GetEntries()) {
    ClonedMCChain = &SignalMCChain;
  } else {
    ClonedMCChain = SignalMCChain.CloneTree(m_Settings.getI("Events_in_MC"));
  }
  RooDataSet MCSignal("MCSignal", "", ClonedMCChain, RooArgList(*m_SignalMBC));
  auto SignalShape = Unique::create<RooKeysPdf*>("SignalShape", "", *m_SignalMBC, MCSignal);
  m_SignalShapeConv = Unique::create<RooFFTConvPdf*>("SignalShapeConv", "", *m_SignalMBC, *SignalShape, *Resolution);
}

void BinnedFitModel::InitializeCombinatorialShape() {
  if(m_Settings.getB("FullyReconstructed")) {
    m_Parameters.insert({"End", Unique::create<RooRealVar*>("End", "", 1.8865)});
    m_Parameters.insert({"c", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_c")});
    m_Combinatorial = Unique::create<RooArgusBG*>("Combinatorial", "", *m_SignalMBC, *m_Parameters["End"], *m_Parameters["c"]);
  } else {
    Chebychev_Shape CombinatorialShape("Combinatorial", m_Settings["MBC_Shape"], m_SignalMBC);
    m_Combinatorial = CombinatorialShape.GetPDF();
  }
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
    } else if(PeakingShape == "CrystalBall") {
      PeakingPDF = new CrystalBall_Shape(Name, m_Settings["MBC_Shape"], m_SignalMBC);
    } else {
      throw std::invalid_argument("Unknown peaking background shape: " + PeakingShape);
    }
    m_PeakingBackgroundShapes.insert({Name, PeakingPDF});
  }
}

RooAddPdf* BinnedFitModel::CreateBinPDF(const std::string &CategoryString) {
  RooArgList Shapes(*m_SignalShapeConv, *m_Combinatorial);
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

void BinnedFitModel::SmearPeakingBackgrounds() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(const auto &Category : m_Category.GetCategories()) {
    for(int i = 0; i < PeakingBackgrounds; i++) {
      std::string Name = Mode + "_PeakingBackground" + std::to_string(i) + "_" + Category;
      if(!m_Settings["MBC_Shape"].contains(Name + "_Yield")) {
	auto YieldVar = static_cast<RooFormulaVar*>(m_Yields[Category + "_PeakingBackground" + std::to_string(i) + "Yield"]);
	double BackgroundSignalRatio = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio");
	double BackgroundSignalRatio_err = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio_err");
	auto BkgSigRatioVar = static_cast<RooRealVar*>(YieldVar->getParameter((Name + "_BackgroundToSignalRatio").c_str()));
	BackgroundSignalRatio += gRandom->Gaus(0.0, BackgroundSignalRatio_err);
	if(BackgroundSignalRatio < 0.0) {
	  BackgroundSignalRatio = 0.0;
	}
	BkgSigRatioVar->setVal(BackgroundSignalRatio);
	if(m_Settings["MBC_Shape"].contains(Name + "_QuantumCorrelationFactor")) {
	  double QCFactor = m_Settings["MBC_Shape"].getD(Name + "_QuantumCorrelationFactor");
	  double QCFactor_err = m_Settings["MBC_Shape"].getD(Name + "_QuantumCorrelationFactor_err");
	  QCFactor += gRandom->Gaus(0.0, QCFactor_err);
	  if(QCFactor < 0.0) {
	    QCFactor = 0.0;
	  }
	  auto QCFactorVar = static_cast<RooRealVar*>(YieldVar->getParameter((Name + "_QuantumCorrelationFactor").c_str()));
	  QCFactorVar->setVal(QCFactor);
	}
      } else {
	double BackgroundYield = m_Settings["MBC_Shape"].getD(Name + "_Yield");
	double BackgroundYield_err = m_Settings["MBC_Shape"].getD(Name + "_Yield_err");
	BackgroundYield += gRandom->Gaus(0.0, BackgroundYield_err);
	if(BackgroundYield < 0.0) {
	  BackgroundYield = 0.0;
	}
	auto PeakBkgYieldVar = static_cast<RooRealVar*>(m_Yields[Category + "_PeakingBackground" + std::to_string(i) + "Yield"]);
	PeakBkgYieldVar->setVal(BackgroundYield);
      }
    }
  }
}
