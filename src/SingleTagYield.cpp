// Martin Duy Tat 4th April 2021

#include<string>
#include<fstream>
#include"TTree.h"
#include"TCanvas.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooArgList.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooGaussian.h"
#include"RooFFTConvPdf.h"
#include"RooAddPdf.h"
#include"RooPlot.h"
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
  m_MBC.setBins(10000, "fft");
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
  //Model.fitTo(Data, PrintEvalErrors(-1));
  Model.fitTo(Data);
  RooPlot *Frame = m_MBC.frame();
  Frame->SetTitle((TagMode + std::string(" Single Tag M_{BC}; M_{BC} (GeV); Events")).c_str());
  Data.plotOn(Frame, Binning(100));
  Model.plotOn(Frame, Color(kBlue));
  Model.plotOn(Frame, Color(kBlue), Components(Argus), LineStyle(kDashed));
  Model.plotOn(Frame, Color(kRed), Components(SignalShapeConv));
  TCanvas c1("c1", "c1", 1600, 1200);
  Frame->Draw();
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
