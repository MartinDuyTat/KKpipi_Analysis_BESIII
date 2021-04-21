// Martin Duy Tat 4th April 2021

#include<string>
#include<fstream>
#include<sstream>
#include"TTree.h"
#include"TCanvas.h"
#include"TPad.h"
#include"TAxis.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooArgList.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooGaussian.h"
#include"RooFFTConvPdf.h"
#include"RooAddPdf.h"
#include"RooExtendPdf.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"SingleTagYield.h"

SingleTagYield::SingleTagYield(TTree *DataTree, TTree *MCSignalTree):
                               m_DataTree(DataTree),
			       m_MCSignalTree(MCSignalTree),
			       m_MBC(RooRealVar("MBC", "MBC", 1.83, 1.885)),
			       m_c(RooRealVar("c", "c", -10, -100, 100)),
			       m_End(RooRealVar("End", "End", 1.8865)),
			       m_Mean(RooRealVar("Mean", "Mean", 0.0, -0.003, 0.003)),
			       m_Sigma(RooRealVar("Sigma", "Sigma", 0.003, 0.00001, 0.010)),
                               m_LuminosityWeight("LuminosityWeight", "LuminosityWeight", 1.0, 0.0, 10.0) {
  m_DataTree->SetBranchStatus("*", 0);
  m_DataTree->SetBranchStatus("MBC", 1);
  m_DataTree->SetBranchStatus("LuminosityWeight", 1);
  m_MCSignalTree->SetBranchStatus("*", 0);
  m_MCSignalTree->SetBranchStatus("MBC", 1);
  m_MBC.setBins(10000, "cache");
  m_MBC.setRange("SignalRange", 1.86, 1.87);
  // Count the number of weighted events
  double N = 0;
  double LuminosityWeight;
  m_DataTree->SetBranchAddress("LuminosityWeight", &LuminosityWeight);
  for(int i = 0; i < m_DataTree->GetEntries(); i++) {
    m_DataTree->GetEntry(i);
    N += LuminosityWeight;
  }
  m_Nsig = RooRealVar("Nsig", "Nsig", N/2, 0.0, N);
  m_Nbkg = RooRealVar("Nbkg", "Nbkg", N/2, 0.0, N);
}

SingleTagYield::~SingleTagYield() {
  for(auto p : m_PeakingMean) {
    delete p;
  }
  for(auto p : m_PeakingSigma) {
    delete p;
  }
  for(auto p : m_PeakingYield) {
    delete p;
  }
  for(auto p : m_PeakingPDF) {
    delete p;
  }
  for(auto p : m_PeakingExPDF) {
    delete p;
  }
}

void SingleTagYield::AddPeakingComponent(const std::string &Filename) {
  std::ifstream Infile(Filename);
  std::string line;
  while(std::getline(Infile, line)) {
    std::stringstream ss(line);
    std::string Name;
    double Mean, Sigma, Yield;
    ss >> Name >> Mean >> Sigma >> Yield;
    m_PeakingName.push_back(Name);
    m_PeakingMean.push_back(new RooRealVar((Name + "Mean").c_str(), (Name + "Mean").c_str(), Mean));
    m_PeakingSigma.push_back(new RooRealVar((Name + "Sigma").c_str(), (Name + "Sigma").c_str(), Sigma));
    m_PeakingYield.push_back(new RooRealVar((Name + "Yield").c_str(), (Name + "Yield").c_str(), Yield));
    m_PeakingPDF.push_back(new RooGaussian(Name.c_str(), Name.c_str(), m_MBC, *m_PeakingMean.back(), *m_PeakingSigma.back()));
    m_PeakingExPDF.push_back(new RooExtendPdf(("Ex" + Name).c_str(), ("Ex" + Name).c_str(), *m_PeakingPDF.back(), *m_PeakingYield.back()));
  }
  Infile.close();
}

void SingleTagYield::FitYield(const std::string &TagMode, const std::string &Filename) {
  using namespace RooFit;
  RooDataSet MCSignal("MCSignal", "MCSignal", m_MCSignalTree, RooArgList(m_MBC));
  RooKeysPdf SignalShape("SignalShape", "SignalShape", m_MBC, MCSignal);
  RooGaussian Resolution("Resolution", "Resolution", m_MBC, m_Mean, m_Sigma);
  RooFFTConvPdf SignalShapeConv("SignalShapeConv", "SignalShapeConv", m_MBC, SignalShape, Resolution);
  RooArgusBG Argus("Argus", "Argus", m_MBC, m_End, m_c);
  RooExtendPdf ExSignalShapeConv("ExSignalShapeConv", "ExSignalShapeConv", SignalShapeConv, m_Nsig, "SignalRange");
  RooExtendPdf ExArgus("ExArgus", "ExArgus", Argus, m_Nbkg);
  RooArgList PDFList(ExSignalShapeConv, ExArgus);
  for(unsigned int i = 0; i < m_PeakingPDF.size(); i++) {
    PDFList.add(*m_PeakingExPDF[i]);
  }
  RooAddPdf Model("Model", "Model", PDFList);
  RooDataSet Data("Data", "Data", m_DataTree, RooArgList(m_MBC, m_LuminosityWeight), "", "LuminosityWeight");
  Model.fitTo(Data, PrintEvalErrors(-1));
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
  Frame->SetTitle((TagMode + std::string(" Single Tag M_{BC}; M_{BC} (GeV); Events")).c_str());
  Data.plotOn(Frame, Binning(100));
  Model.plotOn(Frame, LineColor(kBlue));
  RooHist *Pull = Frame->pullHist();
  Model.plotOn(Frame, LineColor(kBlue), Components(Argus), LineStyle(kDashed));
  Model.plotOn(Frame, LineColor(kRed), Components(SignalShapeConv));
  for(auto iter = m_PeakingName.begin(); iter != m_PeakingName.end(); iter++) {
    Model.plotOn(Frame, LineColor(kGreen), Components(iter->c_str()));
  }
  Frame->Draw();
  Pad2.cd();
  RooPlot *PullFrame = m_MBC.frame();
  PullFrame->addObject(Pull);
  PullFrame->SetMinimum(-5);
  PullFrame->SetMaximum(5);
  PullFrame->SetTitle(";;");
  PullFrame->GetXaxis()->SetLabelFont(0);
  PullFrame->GetXaxis()->SetLabelSize(0);
  PullFrame->GetYaxis()->SetLabelFont(62);
  PullFrame->GetYaxis()->SetLabelSize(0.1);
  PullFrame->Draw();
  c1.cd();
  c1.SaveAs(Filename.c_str());
}

void SingleTagYield::SaveFitParameters(const std::string &Filename) const {
  std::ofstream OutputFile(Filename);
  OutputFile << "Nsig " << m_Nsig.getValV() << " " << m_Nsig.getError() << "\n";
  OutputFile << "Nbkg " << m_Nbkg.getValV() << " " << m_Nbkg.getError() << "\n";
  OutputFile << "Sigma " << m_Sigma.getValV() << " " << m_Sigma.getError() << "\n";
  OutputFile << "Mean " << m_Mean.getValV() << " " << m_Mean.getError() << "\n";
  OutputFile << "c " << m_c.getValV() << " " << m_c.getError() << "\n";
  OutputFile.close();
}

double SingleTagYield::GetYield() const {
  return m_Nsig.getVal();
}
