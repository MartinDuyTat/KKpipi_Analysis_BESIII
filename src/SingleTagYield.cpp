// Martin Duy Tat 4th April 2021

#include<string>
#include<utility>
#include<fstream>
#include<stdexcept>
#include<numeric>
#include<memory>
#include<vector>
#include"TTree.h"
#include"TCanvas.h"
#include"TPad.h"
#include"TAxis.h"
#include"TLine.h"
#include"TH1D.h"
#include"TFile.h"
#include"TRandom.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooDataHist.h"
#include"RooArgList.h"
#include"RooArgSet.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooGaussian.h"
#include"RooFFTConvPdf.h"
#include"RooAddPdf.h"
#include"RooFormulaVar.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"RooStats/SPlot.h"
#include"RooStats/RooStatsUtils.h"
#include"SingleTagYield.h"
#include"Settings.h"
#include"Unique.h"
#include"Utilities.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"
#include"RooShapes/CrystalBall_Shape.h"

SingleTagYield::SingleTagYield(TTree *DataTree, TTree *MCSignalTree, const Settings &settings):
                               m_DataTree(DataTree),
			       m_MCSignalTree(MCSignalTree),
			       m_Settings(settings),
			       m_MBC("MBC", "", 1.83, 1.8865),
                               m_LuminosityWeight("LuminosityWeight", "", 1.0, 0.0, 10.0) {
  m_DataTree->SetBranchStatus("*", 0);
  m_DataTree->SetBranchStatus("MBC", 1);
  m_DataTree->SetBranchStatus("LuminosityWeight", 1);
  m_MBC.setBins(1000, "cache");
  InitializeSignalShape();
  InitializeArgus();
  InitializePeakingBackgrounds();
  InitializeFitShape();
}

SingleTagYield::~SingleTagYield() {
  for(auto &PeakingBackground : m_PeakingBackgrounds) {
    delete PeakingBackground;
  }
}

void SingleTagYield::InitializeSignalShape() {
  std::string Name = m_Settings.get("Mode") + "_SingleTag_";
  RooAbsPdf* Resolution = nullptr;
  if(m_Settings["MBC_Shape"].contains(Name + "SingleGaussianResolution") && 
     m_Settings["MBC_Shape"].getB(Name + "SingleGaussianResolution")) {
    auto Mean = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Mean");
    m_Parameters.insert({"Mean", Mean});
    auto Sigma = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Sigma");
    m_Parameters.insert({"Sigma", Sigma});
    Resolution = Unique::create<RooGaussian*>("Gaussian", "", m_MBC, *m_Parameters["Mean"], *m_Parameters["Sigma"]);
  } else {
    m_Parameters.insert({"frac", Unique::create<RooRealVar*>((Name + "frac").c_str(), "", 0.5, 0.0, 1.0)});
    auto Mean1 = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Mean1");
    auto Mean2 = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Mean2");
    m_Parameters.insert({"Mean1", Mean1});
    m_Parameters.insert({"Mean2", Mean2});
    auto Sigma1 = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Sigma1");
    auto Sigma2 = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Sigma2");
    m_Parameters.insert({"Sigma1", Sigma1});
    m_Parameters.insert({"Sigma2", Sigma2});
    auto Gaussian1 = Unique::create<RooGaussian*>("Gaussian1", "", m_MBC, *m_Parameters["Mean1"], *m_Parameters["Sigma1"]);
    auto Gaussian2 = Unique::create<RooGaussian*>("Gaussian2", "", m_MBC, *m_Parameters["Mean2"], *m_Parameters["Sigma2"]);
    Resolution = Unique::create<RooAddPdf*>("Resolution", "", RooArgList(*Gaussian1, *Gaussian2), *m_Parameters["frac"]);
  }
  RooDataSet MCSignal("MCSignal", "", m_MCSignalTree, RooArgList(m_MBC));
  auto SignalShape = Unique::create<RooKeysPdf*>("SignalShape", "", m_MBC, MCSignal);
  auto SignalShapeConv = Unique::create<RooFFTConvPdf*>("SignalShapeConv", "", m_MBC, *SignalShape, *Resolution);
  m_ModelPDFs.add(*SignalShapeConv);
  auto SingleTag_Yield = Utilities::load_param(m_Settings["MBC_Shape"], Name + "Yield");
  m_ModelYields.add(*SingleTag_Yield);
  m_Parameters.insert({"Yield", SingleTag_Yield});
}

void SingleTagYield::InitializeArgus() {
  m_Parameters.insert({"End", Unique::create<RooRealVar*>("End", "", 1.8865)});
  auto c = Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_SingleTag_c");
  m_Parameters.insert({"c", c});
  auto Argus = Unique::create<RooArgusBG*>("Argus", "", m_MBC, *m_Parameters["End"], *m_Parameters["c"]);
  m_ModelPDFs.add(*Argus);
  auto SingleTag_BackgroundYield = Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_SingleTag_CombinatorialYield");
  m_ModelYields.add(*SingleTag_BackgroundYield);
  m_Parameters.insert({"CombinatorialYield", SingleTag_BackgroundYield});
}

void SingleTagYield::InitializePeakingBackgrounds() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    std::string PeakingShape = m_Settings["MBC_Shape"].get(Name + "_Shape");
    FitShape *PeakingPDF = nullptr;
    if(PeakingShape == "DoubleGaussian") {
      PeakingPDF = new DoubleGaussian_Shape(Name, m_Settings["MBC_Shape"], &m_MBC);
    } else if(PeakingShape == "DoubleCrystalBall") {
      PeakingPDF = new DoubleCrystalBall_Shape(Name, m_Settings["MBC_Shape"], &m_MBC);
    } else if(PeakingShape == "CrystalBall") {
      PeakingPDF = new CrystalBall_Shape(Name, m_Settings["MBC_Shape"], &m_MBC);
    } else {
      throw std::invalid_argument("Unknown peaking background shape");
    }
    if(m_Settings["MBC_Shape"].contains(Name + "_yield")) {
      double PeakingYield = m_Settings["MBC_Shape"].getD(Name + "_yield");
      std::cout << "Adding peaking background with yield: " << PeakingYield << "\n";
      // This is automatically handled in FitShape
    } else {
      double BackgroundToSignalRatio = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio");
      PeakingPDF->UseRelativeYield(m_Parameters["Yield"], BackgroundToSignalRatio);
      std::cout << "Adding peaking background with background-to-signal ratio: " << BackgroundToSignalRatio << "\n";
    }
    m_ModelPDFs.add(*PeakingPDF->GetPDF());
    m_ModelYields.add(*PeakingPDF->GetYield());
    m_PeakingBackgrounds.push_back(PeakingPDF);
  }
}

void SingleTagYield::InitializeFitShape() {
  m_FullModel = Unique::create<RooAddPdf*>("MBC_Model", "", m_ModelPDFs, m_ModelYields);
} 

void SingleTagYield::FitYield() {
  using namespace RooFit;
  RooArgList Variables(m_MBC, m_LuminosityWeight);
  std::string MassCut(""), MassCutBinned("");
  std::unique_ptr<RooRealVar> InvMassVar;
  if(m_Settings.contains("InvariantMassVariable")) {
    std::string MassVarName = m_Settings.get("InvariantMassVariable");
    m_DataTree->SetBranchStatus(MassVarName.c_str(), 1);
    double LowMassCut = m_Settings.getD("InvariantMassVariable_low");
    double HighMassCut = m_Settings.getD("InvariantMassVariable_high");
    InvMassVar = std::unique_ptr<RooRealVar>(new RooRealVar(MassVarName.c_str(), "", LowMassCut, HighMassCut));
    MassCut = MassVarName + " > " + std::to_string(LowMassCut) + " && " + MassVarName + " < " + std::to_string(HighMassCut);
    MassCutBinned = MassCut + "*LuminosityWeight";
    Variables.add(*InvMassVar);
  }
  TH1D h1("h1", "h1", m_Settings.getI("Bins_in_fit"), 1.83, 1.8865);
  m_DataTree->Draw("MBC >> h1", MassCutBinned.c_str(), "goff");
  RooDataHist BinnedData("BinnedData", "BinnedData", RooArgList(m_MBC), &h1);
  RooDataSet Data("Data", "Data", m_DataTree, Variables, MassCut.c_str(), "LuminosityWeight");
  RooArgSet *Parameters = m_FullModel->getParameters(m_MBC);
  m_InitialParameters = Parameters->snapshot();
  if(m_Settings.get("FitType") != "NoFit") {
    m_Result = m_FullModel->fitTo(BinnedData, Save(), Strategy(2));
    if(m_Settings.get("FitType") == "UnbinnedFit") {
      m_Result = m_FullModel->fitTo(Data, Save(), Strategy(2), NumCPU(4));
    }
    SaveFitParameters();
  }
  PlotSingleTagYield(Data);
  if(m_Settings.getB("YieldSystematics")) {
    double SystError;
    int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(m_Settings.get("Mode") + "_PeakingBackgrounds");
    if(PeakingBackgrounds > 0) {
      std::vector<double> FittedYields;
      gRandom->SetSeed(m_Settings.getI("Seed"));
      for(int i = 0; i < m_Settings.getI("NumberRuns"); i++) {
	std::cout << "Starting systematics fit number: " << i << "\n";
	*Parameters = *m_InitialParameters;
	SmearPeakingBackgrounds();
	auto Result = m_FullModel->fitTo(BinnedData, Strategy(2), Save());
	Result->Print("V");
	FittedYields.push_back(m_Parameters["Yield"]->getVal());
      }
      SystError = TMath::RMS(FittedYields.begin(), FittedYields.end());
    } else {
      SystError = 0.0;
    }
    std::ofstream OutputFile(m_Settings.get("ResultsFilename"), std::ios_base::app);
    std::string YieldName = m_Settings.get("Mode") + "_SingleTag_Yield_PeakingBackgrounds";
    OutputFile << YieldName << "_syst_err " << SystError << "\n";
    OutputFile.close();
  }
  if(m_Settings.contains("sPlotReweight") && m_Settings.getB("sPlotReweight")) {
    sPlotReweight(Data);
  }
}

void SingleTagYield::PlotSingleTagYield(const RooDataSet &Data) const {
  using namespace RooFit;
  TCanvas c1("c1", "c1", 1600, 1200);
  TPad Pad1("Pad1", "Pad1", 0.0, 0.25, 1.0, 1.0);
  TPad Pad2("Pad2", "Pad2", 0.0, 0.0, 1.0, 0.25);
  Pad1.Draw();
  Pad2.Draw();
  Pad1.SetBottomMargin(0.1);
  Pad1.SetTopMargin(0.1);
  Pad1.SetBorderMode(0);
  Pad2.SetBorderMode(0);
  Pad2.SetBottomMargin(0.1);
  Pad2.SetTopMargin(0.05);
  Pad1.cd();
  RooPlot *Frame = m_MBC.frame();
  std::string TagMode = m_Settings.get("Mode");
  Frame->SetTitle((TagMode + std::string(" Single Tag M_{BC}; M_{BC} (GeV); Events")).c_str());
  Data.plotOn(Frame, Binning(200));
  double TotalEvents = std::accumulate(m_ModelYields.begin(), m_ModelYields.end(), 0.0, [] (double a, const auto &b) { return a + static_cast<RooRealVar*>(b)->getVal(); });
  m_FullModel->plotOn(Frame, LineColor(kBlue), Normalization(TotalEvents, RooAbsReal::NumEvent));
  RooHist *Pull = Frame->pullHist();
  m_FullModel->plotOn(Frame, LineColor(kBlue), Components("Argus"), LineStyle(kDashed), Normalization(TotalEvents, RooAbsReal::NumEvent));
  m_FullModel->plotOn(Frame, LineColor(kRed), Components("SignalShapeConv"), Normalization(TotalEvents, RooAbsReal::NumEvent));
  if(m_PeakingBackgrounds.size() > 0) {
    std::string PeakingList = m_PeakingBackgrounds[0]->GetPDF()->GetName();
    for(unsigned int i = 1; i < m_PeakingBackgrounds.size(); i++) {
      PeakingList += std::string(",") + m_PeakingBackgrounds[i]->GetPDF()->GetName();
    }
    m_FullModel->plotOn(Frame, LineColor(kMagenta), Components(PeakingList.c_str()));
  }
  Frame->Draw();
  Pad2.cd();
  RooPlot *PullFrame = m_MBC.frame();
  PullFrame->addObject(Pull);
  TLine *Line = new TLine(1.83, 0.0, 1.8865, 0.0);
  PullFrame->addObject(Line);
  PullFrame->SetMinimum(-5);
  PullFrame->SetMaximum(5);
  PullFrame->SetTitle(";;");
  PullFrame->GetXaxis()->SetLabelFont(0);
  PullFrame->GetXaxis()->SetLabelSize(0);
  PullFrame->GetYaxis()->SetLabelFont(62);
  PullFrame->GetYaxis()->SetLabelSize(0.1);
  PullFrame->Draw();
  c1.cd();
  std::string PlotFilename = m_Settings.get("MBCPlotFilename");
  c1.SaveAs(PlotFilename.c_str());
}

void SingleTagYield::SaveFitParameters() const {
  m_Result->Print("V");
  std::string Mode = m_Settings.get("Mode");
  std::ofstream OutputFile(m_Settings.get("ResultsFilename"));
  OutputFile << "status " << m_Result->status() << "\n";
  OutputFile << "covQual " << m_Result->covQual() << "\n";
  for(const auto Param : m_Parameters) {
    OutputFile << Mode + "_" + Param.first << " " << Param.second->getVal() << "\n";
    OutputFile << Mode + "_" + Param.first << "_err " << Param.second->getError() << "\n";
  }
  auto Yield = m_Parameters.at("Yield");
  OutputFile << Mode + "_SingleTag_Yield " << Yield->getVal() << "\n";
  OutputFile << Mode + "_SingleTag_Yield_err " << Yield->getError() << "\n";
  OutputFile.close();
}

std::pair<double, double> SingleTagYield::CalculateSingleTagYield() const {
  using namespace RooFit;
  return std::pair<double, double>{m_Parameters.at("Yield")->getVal(), m_Parameters.at("Yield")->getPropagatedError(*m_Result)};
}

void SingleTagYield::sPlotReweight(RooDataSet &Data) {
  auto FloatingParameters = m_FullModel->getParameters(Data);
  for(auto Parameter : *FloatingParameters) {
    static_cast<RooRealVar*>(Parameter)->setConstant(true);
  }
  RooArgList sPlotYields;
  for(auto YieldParameter : m_ModelYields) {
    auto YieldParameter_Var = dynamic_cast<RooRealVar*>(YieldParameter);
    if(YieldParameter_Var) {
      YieldParameter_Var->setConstant(false);
      sPlotYields.add(*YieldParameter_Var);
    }
  }
  RooStats::SPlot("sData", "", Data, m_FullModel, sPlotYields);
  TFile Outfile(m_Settings.get("sPlotFilename").c_str(), "RECREATE");
  auto Tree = RooStats::GetAsTTree(m_Settings.get("TreeName").c_str(), m_Settings.get("TreeName").c_str(), Data);
  Tree->Write();
  Outfile.Close();
}

void SingleTagYield::SmearArgusEndPoint() {
  m_Parameters["End"]->setVal(1.8865 + gRandom->Gaus(0.0, 0.0005));
}

void SingleTagYield::SmearPeakingBackgrounds() {
  std::string Mode = m_Settings.get("Mode");
  int PeakingBackgrounds = m_Settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    auto YieldVar = m_PeakingBackgrounds[i]->GetYield();
    if(m_Settings["MBC_Shape"].contains(Name + "_yield")) {
      double PeakingYield = m_Settings["MBC_Shape"].getD(Name + "_yield");
      double PeakingYield_err = m_Settings["MBC_Shape"].getD(Name + "_yield_err");
      PeakingYield += gRandom->Gaus(0.0, PeakingYield_err);
      if(PeakingYield <= 0.0) {
	PeakingYield = 0.0;
      }
      static_cast<RooRealVar*>(YieldVar)->setVal(PeakingYield);
    } else {
      double BackgroundToSignalRatio = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio");
      double BackgroundToSignalRatio_err = m_Settings["MBC_Shape"].getD(Name + "_BackgroundToSignalRatio_err");
      BackgroundToSignalRatio += gRandom->Gaus(0.0, BackgroundToSignalRatio_err);
      if(BackgroundToSignalRatio <= 0.0) {
	BackgroundToSignalRatio = 0.0;
      }
      auto BkgToSigVar = static_cast<RooFormulaVar*>(YieldVar)->getParameter((Name + "_BackgroundToSignalYieldRatio").c_str());
      static_cast<RooRealVar*>(BkgToSigVar)->setVal(BackgroundToSignalRatio);
    }
  }
}
