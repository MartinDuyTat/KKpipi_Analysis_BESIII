// Martin Duy Tat 4th April 2021

#include<string>
#include<fstream>
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
#include"RooPlot.h"
#include"RooHist.h"
#include"SingleTagYield.h"

SingleTagYield::SingleTagYield(TTree *DataTree, TTree *MCSignalTree):
                               m_DataTree(DataTree),
			       m_MCSignalTree(MCSignalTree),
			       m_MBC(RooRealVar("MBC", "MBC", 1.83, 1.89)),
			       m_c(RooRealVar("c", "c", -10, -100, 100)),
			       m_End(RooRealVar("End", "End", 1.8865)),
			       m_Mean(RooRealVar("Mean", "Mean", 0.0)),
			       m_Sigma(RooRealVar("Sigma", "Sigma", 0.003, 0.0001, 0.100)),
			       m_Nsig(RooRealVar("Nsig", "Nsig", m_DataTree->GetEntries()/2, 0.0, m_DataTree->GetEntries())),
			       m_Nbkg(RooRealVar("Nbkg", "Nbkg", m_DataTree->GetEntries()/2, 0.0, m_DataTree->GetEntries())) {
  m_DataTree->SetBranchStatus("*", 0);
  m_DataTree->SetBranchStatus("MBC", 1);
  m_MCSignalTree->SetBranchStatus("*", 0);
  m_MCSignalTree->SetBranchStatus("MBC", 1);
  m_MBC.setBins(10000, "cache");
}

void SingleTagYield::FitYield(const std::string &TagMode, const std::string &Filename) {
  using namespace RooFit;
  RooDataSet MCSignal("MCSignal", "MCSignal", m_MCSignalTree, RooArgList(m_MBC));
  RooKeysPdf SignalShape("SignalShape", "SignalShape", m_MBC, MCSignal);
  RooGaussian Resolution("Resolution", "Resolution", m_MBC, m_Mean, m_Sigma);
  RooFFTConvPdf SignalShapeConv("SignalShapeConv", "SignalShapeConv", m_MBC, SignalShape, Resolution);
  RooArgusBG Argus("Argus", "Argus", m_MBC, m_End, m_c);
  RooAddPdf Model("Model", "Model", RooArgList(SignalShapeConv, Argus), RooArgList(m_Nsig, m_Nbkg));
  RooDataSet Data("Data", "Data", m_DataTree, RooArgList(m_MBC));
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
  OutputFile << "c " << m_c.getValV() << " " << m_c.getError() << "\n";
  OutputFile.close();
}

double SingleTagYield::GetYield() const {
  return m_Nsig.getVal();
}
