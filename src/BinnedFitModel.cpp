// Martin Duy Tat 25th November 2021

#include<string>
#include<stdexcept>
#include"TTree.h"
#include"TRandom.h"
#include"TFile.h"
#include"TROOT.h"
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
  std::string PDFName = "Simultaneous_PDF_KKpipi_vs_" + settings.get("Mode");
  m_PDF = new RooSimultaneous(PDFName.c_str(), "", *CategoryVariable);
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
  m_SignalYields.setName((Mode + "_SignalYields").c_str());
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(const auto &CategoryString : m_Category.GetCategories()) {
    std::string CombinatorialName = CategoryString + "_CombinatorialYield";
    m_Yields.insert({
      CombinatorialName,
      Unique::create<RooRealVar*>(CombinatorialName.c_str(), "",
				  1.0, 0.0, 1000.0)});
    std::string SignalName = CategoryString + "_SignalYield";
    m_Yields.insert({
      SignalName,
      Unique::create<RooRealVar*>(SignalName.c_str(), "",
				  10.0, 0.0, 1000.0)});
    m_SignalYields.add(*m_Yields.at(SignalName));
    for(int i = 0; i < PeakingBackgrounds; i++) {
      std::string Name = Mode + "_PeakingBackground";
      Name += std::to_string(i) + "_" + CategoryString;
      if(!m_Settings["MBC_Shape"].contains(Name + "_Yield")) {
	auto BackgroundSignalRatio =
	  Utilities::load_param(m_Settings["BackgroundToSignalRatio"],
				Name + "_BackgroundToSignalRatio");
	std::cout << CategoryString << " peaking background " << i;
	std::cout << " background-to-signal ratio: ";
	std::cout << BackgroundSignalRatio->getVal() << "\n";
	// Quantum correlation is accounted for with a correction factor
	RooRealVar *QCFactor = nullptr;
	if(m_Settings.contains_subsettings("QuantumCorrelationFactor")) {
	  QCFactor = Utilities::load_param(m_Settings["QuantumCorrelationFactor"],
					   Name + "_QuantumCorrelationFactor");
	  std::cout << CategoryString << " peaking background " << i;
	  std::cout << " quantum correlation correction factor: ";
	  std::cout << QCFactor->getVal() << "\n";
	} else {
	  QCFactor =
	    Unique::create<RooRealVar*>((Name + "QuantumCorrelationFactor").c_str(),
					"", 1.0);
	  std::cout << CategoryString << " peaking background " << i;
	  std::cout << " quantum correlation correction factor: None\n";
	}
	RooArgList PeakingYieldParameters(*m_Yields.at(SignalName),
					  *BackgroundSignalRatio,
					  *QCFactor);
	RooAbsReal *PeakingYield =
	  Unique::create<RooFormulaVar*>((Name + "_Yield").c_str(),
					 "@0*@1*@2",
					 PeakingYieldParameters);
	std::string PeakingYieldName = CategoryString + "_PeakingBackground";
	PeakingYieldName += std::to_string(i) + "Yield";
	m_Yields.insert({PeakingYieldName, PeakingYield});
      } else {
	double BackgroundYield = m_Settings["MBC_Shape"].getD(Name + "_Yield");
	RooAbsReal *PeakingYield =
	  Unique::create<RooRealVar*>((Name + "_Yield").c_str(),
				      "", BackgroundYield);
	std::cout << CategoryString << " peaking background " << i << " yield: ";
	std::cout << BackgroundYield << "\n";
	std::string PeakingYieldName = CategoryString + "_PeakingBackground";
	PeakingYieldName += std::to_string(i) + "Yield";
	m_Yields.insert({PeakingYieldName, PeakingYield});
      }
    }
  }
}

void BinnedFitModel::InitializeSignalShape() {
  m_Parameters.insert({
    "Mean1",
    Utilities::load_param(m_Settings["MBC_Shape"],
			  m_Settings.get("Mode") + "_DoubleTag_Mean1")});
  m_Parameters.insert({
    "Sigma1",
    Utilities::load_param(m_Settings["MBC_Shape"],
			  m_Settings.get("Mode") + "_DoubleTag_Sigma1")});
  auto Resolution = Unique::create<RooGaussian*>("Gaussian1",
						 "",
						 *m_SignalMBC,
						 *m_Parameters["Mean1"],
						 *m_Parameters["Sigma1"]);
  TChain SignalMCChain(m_Settings.get("TreeName").c_str());
  std::string SignalMCFilename =
    m_Settings["Datasets_WithDeltaECuts"].get("SignalMC_DT");
  SignalMCFilename = Utilities::ReplaceString(SignalMCFilename,
					      "TAG",
					      m_Settings.get("Mode"));
  SignalMCChain.Add(SignalMCFilename.c_str());
  SignalMCChain.SetBranchStatus("*", 0);
  SignalMCChain.SetBranchStatus(m_Settings.get("FitVariable").c_str(), 1);
  TTree *ClonedMCChain = nullptr;
  gROOT->cd();
  if(m_Settings.getI("Events_in_MC") < 0 ||
     m_Settings.getI("Events_in_MC") > SignalMCChain.GetEntries()) {
    ClonedMCChain = SignalMCChain.CloneTree();
  } else {
    ClonedMCChain = SignalMCChain.CloneTree(m_Settings.getI("Events_in_MC"));
  }
  RooDataSet MCSignal("MCSignal", "", ClonedMCChain, RooArgList(*m_SignalMBC));
  auto SignalShape = Unique::create<RooKeysPdf*>("SignalShape",
						 "",
						 *m_SignalMBC,
						 MCSignal);
  m_SignalShapeConv = Unique::create<RooFFTConvPdf*>("SignalShapeConv",
						     "",
						     *m_SignalMBC,
						     *SignalShape,
						     *Resolution);
}

void BinnedFitModel::InitializeCombinatorialShape() {
  if(m_Settings.getB("FullyReconstructed")) {
    m_Parameters.insert({"End", Unique::create<RooRealVar*>("End", "", 1.8865)});
    m_Parameters.insert({
      "c",
      Utilities::load_param(m_Settings["MBC_Shape"],
			    m_Settings.get("Mode") + "_DoubleTag_c")});
    m_Combinatorial = Unique::create<RooArgusBG*>("Combinatorial",
						  "",
						  *m_SignalMBC,
						  *m_Parameters["End"],
						  *m_Parameters["c"]);
  } else {
    Chebychev_Shape CombinatorialShape("Combinatorial",
				       m_Settings["MBC_Shape"],
				       m_SignalMBC);
    m_Combinatorial = CombinatorialShape.GetPDF();
  }
}

void BinnedFitModel::InitializePeakingBackgroundShapes() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    std::string PeakingShape = m_Settings["MBC_Shape"].get(Name + "_Shape");
    FitShape *PeakingPDF = nullptr;
    if(PeakingShape == "DoubleGaussian") {
      PeakingPDF = new DoubleGaussian_Shape(Name,
					    m_Settings["MBC_Shape"],
					    m_SignalMBC);
    } else if(PeakingShape == "DoubleCrystalBall") {
      PeakingPDF = new DoubleCrystalBall_Shape(Name,
					       m_Settings["MBC_Shape"],
					       m_SignalMBC);
    } else if(PeakingShape == "CrystalBall") {
      PeakingPDF = new CrystalBall_Shape(Name,
					 m_Settings["MBC_Shape"],
					 m_SignalMBC);
    } else {
      std::string Invalid = "Unknown peaking background shape: " + PeakingShape;
      throw std::invalid_argument(Invalid);
    }
    m_PeakingBackgroundShapes.insert({Name, PeakingPDF});
  }
}

RooAddPdf* BinnedFitModel::CreateBinPDF(const std::string &CategoryString) {
  RooArgList Shapes(*m_SignalShapeConv, *m_Combinatorial);
  RooArgList Yields(*m_Yields.at(CategoryString + "_SignalYield"),
		    *m_Yields.at(CategoryString + "_CombinatorialYield"));
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    Shapes.add(*m_PeakingBackgroundShapes.at(Name)->GetPDF());
    std::string PeakingYieldName = CategoryString + "_PeakingBackground";
    PeakingYieldName += std::to_string(i) + "Yield";
    Yields.add(*m_Yields.at(PeakingYieldName));
  }
  return Unique::create<RooAddPdf*>((CategoryString + "_PDF").c_str(),
				    "",
				    Shapes,
				    Yields);
}

void BinnedFitModel::InitializePDF() {
  for(const auto &CategoryString : m_Category.GetCategories()) {
    RooAddPdf* BinPDF = CreateBinPDF(CategoryString);
    m_PDF->addPdf(*BinPDF, CategoryString.c_str());
  }
}

void BinnedFitModel::SetGeneratorYields() {
  for(const auto &CategoryString : m_Category.GetCategories()) {
    std::string SignalName = CategoryString + "_SignalYield";
    double GeneratorYield = m_Settings["GeneratorYields"].getD(CategoryString);
    static_cast<RooRealVar*>(m_Yields[SignalName])->setVal(GeneratorYield);
  }
}

void BinnedFitModel::PrepareSmearing() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name(Mode + "_PeakingBackground" + std::to_string(i));
    if(m_Settings["MBC_Shape"].getB(Name + "_Correlated")) {
      std::cout << Mode << " peaking background " << i;
      std::cout << ": Will use Cholesky decomposition\n";
      std::string Filename = Name + "_BackgroundToSignalRatio_CovMatrix.root";
      TFile BkgSigRatioFile(Filename.c_str(), "READ");
      TMatrixT<double> *BkgSigRatioCovMatrix = nullptr;
      BkgSigRatioFile.GetObject("CovMatrix", BkgSigRatioCovMatrix);
      m_CholeskyDecompositions.insert({Name + "_BackgroundToSignalRatio",
	                               CholeskySmearing(*BkgSigRatioCovMatrix)});
      Filename = Name + "_QuantumCorrelationFactor_CovMatrix.root";
      TFile QCFactorFile(Filename.c_str(), "READ");
      TMatrixT<double> *QCFactorCovMatrix = nullptr;
      QCFactorFile.GetObject("CovMatrix", QCFactorCovMatrix);
      m_CholeskyDecompositions.insert({Name + "_QuantumCorrelationFactor",
                                       CholeskySmearing(*QCFactorCovMatrix)});
    }
  }
}

void BinnedFitModel::SmearPeakingBackgrounds() {
  // First smear the correlated peaking backgrounds, if any
  for(auto &CholeskyDecomposition : m_CholeskyDecompositions) {
    CholeskyDecomposition.second.Smear();
  }
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  // Loop over all peaking backgrounds
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string BackgroundName(Mode + "_PeakingBackground" + std::to_string(i));
    // Check if this background has correlated backgrounds
    std::string BkgToSigName = BackgroundName + "_BackgroundToSignalRatio";
    bool Correlated = m_CholeskyDecompositions.find(BkgToSigName) !=
                      m_CholeskyDecompositions.end();
    // Loop over all bins
    for(const auto &Category : m_Category.GetCategories()) {
      std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
      Name += "_" + Category;
      if(!m_Settings["MBC_Shape"].contains(Name + "_Yield")) {
	// If peaking background is expressed as a background-to-signal ratio
	// with quantum correlation correction
	std::string PeakingYieldName = Category + "_PeakingBackground";
	PeakingYieldName += std::to_string(i) + "Yield";
	auto YieldVar =
	  static_cast<RooFormulaVar*>(m_Yields[PeakingYieldName]);
	double BackgroundSignalRatio =
	  m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio");
	if(Correlated) {
	  // If peaking background is correlated, get smearing from Cholesky decomposition
	  int CategoryIndex = m_Category.GetCategoryIndex(Category);
	  BackgroundSignalRatio +=
	    m_CholeskyDecompositions.at(BkgToSigName).GetSmearing(CategoryIndex);
	} else {
	  // Else just generate smearing
	  double BackgroundSignalRatio_err =
	    m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio_err");
	  BackgroundSignalRatio += gRandom->Gaus(0.0, BackgroundSignalRatio_err);
	}
	if(BackgroundSignalRatio < 0.0) {
	  BackgroundSignalRatio = 0.0;
	}
	// Update RooFit variable with smeared yield
	auto BkgSigRatioParameter =
	  YieldVar->getParameter((Name + "_BackgroundToSignalRatio").c_str());
	auto BkgSigRatioVar = static_cast<RooRealVar*>(BkgSigRatioParameter);
	BkgSigRatioVar->setVal(BackgroundSignalRatio);
	// Repeat the same with the quantum correlation factors
	if(m_Settings["MBC_Shape"].contains(Name + "_QuantumCorrelationFactor")) {
	  double QCFactor =
	    m_Settings["MBC_Shape"].getD(Name + "_QuantumCorrelationFactor");
	  if(Correlated) {
	    std::string BkgQCName = BackgroundName + "_QuantumCorrelationFactor";
	    int CategoryIndex = m_Category.GetCategoryIndex(Category);
	    QCFactor +=
	      m_CholeskyDecompositions.at(BkgQCName).GetSmearing(CategoryIndex);
	  } else {
	    double QCFactor_err =
	      m_Settings["MBC_Shape"].getD(Name + "_QuantumCorrelationFactor_err");
	    QCFactor += gRandom->Gaus(0.0, QCFactor_err);
	  }
	  if(QCFactor < 0.0) {
	    QCFactor = 0.0;
	  }
	  auto QCFactorParameter =
	    YieldVar->getParameter((Name + "_QuantumCorrelationFactor").c_str());
	  auto QCFactorVar = static_cast<RooRealVar*>(QCFactorParameter);
	  QCFactorVar->setVal(QCFactor);
	}
      } else {
	// If a peaking background yield is given
	double BackgroundYield = m_Settings["MBC_Shape"].getD(Name + "_Yield");
	double BackgroundYield_err =
	  m_Settings["MBC_Shape"].getD(Name + "_Yield_err");
	BackgroundYield += gRandom->Gaus(0.0, BackgroundYield_err);
	if(BackgroundYield < 0.0) {
	  BackgroundYield = 0.0;
	}
	auto PeakBkgYieldParameter =
	  m_Yields[Category + "_PeakingBackground" + std::to_string(i) + "Yield"];
	auto PeakBkgYieldVar = static_cast<RooRealVar*>(PeakBkgYieldParameter);
	PeakBkgYieldVar->setVal(BackgroundYield);
      }
    }
  }
}
