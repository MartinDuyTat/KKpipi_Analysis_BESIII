// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include"TTree.h"
#include"TH1D.h"
#include"TCanvas.h"
#include"TPad.h"
#include"RooRealVar.h"
#include"RooDataHist.h"
#include"RooDataSet.h"
#include"RooArgList.h"
#include"RooArgList.h"
#include"RooGaussian.h"
//#include"RooPolynomial.h"
#include"RooGenericPdf.h"
#include"RooAddPdf.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"DeltaEFit.h"

DeltaEFit::DeltaEFit(TTree *Tree):
                     m_Tree(Tree),
		     m_DeltaE(RooRealVar("DeltaE", "DeltaE", -0.08, 0.08)),
		     m_Mean1(RooRealVar("Mean1", "Mean1", 0.0, -0.002, 0.002)),
		     m_Mean2(RooRealVar("Mean2", "Mean2", 0.0, -0.002, 0.002)),
		     m_Sigma1(RooRealVar("Sigma1", "Sigma1", 0.001, 0.00001, 0.02)),
                     m_Sigma2(RooRealVar("Sigma2", "Sigma2", 0.001, 0.00001, 0.02)),
                     m_a1(RooRealVar("a1", "a1", 0.0, -100.0, 100.0)),
                     m_b1(RooRealVar("b1", "b1", 0.0, -1000.0, 1000.0)),
                     m_a2(RooRealVar("a2", "a2", 0.0, -100.0, 100.0)),
                     m_b2(RooRealVar("b2", "b2", 0.0, -1000.0, 1000.0)),
                     m_LuminosityWeight("LuminosityWeight", "LuminosityWeight", 1.0, 0.0, 10.0) {
  m_Tree->SetBranchStatus("*", 0);
  m_Tree->SetBranchStatus("DeltaE" , 1);
  m_Tree->SetBranchStatus("DataSetType", 1);
  m_Tree->SetBranchStatus("IsSameDMother", 1);
  m_Tree->SetBranchStatus("PIDTrue", 1);
  m_Tree->SetBranchStatus("LuminosityWeight", 1);
  // Count the number of weighted events
  double N = 0;
  double LuminosityWeight;
  m_Tree->SetBranchAddress("LuminosityWeight", &LuminosityWeight);
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    N += LuminosityWeight;
  }
  m_Nsig1 = RooRealVar("Nsig1", "Nsig1", N/2, 0.0, N);
  m_Nsig2 = RooRealVar("Nsig2", "Nsig2", N/3, 0.0, N);
  m_Nbkg = RooRealVar("Nbkg", "Nbkg", N/2, 0.0, N);
}

void DeltaEFit::SetParameters(const std::string &Filename) {
  std::ifstream Infile(Filename);
  std::string ParameterName;
  double Value, Minimum, Maximum;
  while(Infile.peek() != EOF) {
    Infile >> ParameterName >> Value >> Minimum >> Maximum;
    if(ParameterName == "Nsig1") {
      m_Nsig1 = RooRealVar("Nsig1", "Nsig1", Value, Minimum, Maximum);
    } else if(ParameterName == "Nsig2") {
      m_Nsig2 = RooRealVar("Nsig2", "Nsig2", Value, Minimum, Maximum);
    } else if(ParameterName == "Nbkg") {
      m_Nbkg = RooRealVar("Nbkg", "Nbkg", Value, Minimum, Maximum);
    } else if(ParameterName == "Mean1") {
      m_Mean1 = RooRealVar("Mean1", "Mean1", Value, Minimum, Maximum);
    } else if(ParameterName == "Mean2") {
      m_Mean2 = RooRealVar("Mean2", "Mean2", Value, Minimum, Maximum);
    } else if(ParameterName == "Sigma1") {
      m_Sigma1 = RooRealVar("Sigma1", "Sigma1", Value, Minimum, Maximum);
    } else if(ParameterName == "Sigma2") {
      m_Sigma2 = RooRealVar("Sigma2", "Sigma2", Value, Minimum, Maximum);
    } else if(ParameterName == "a1") {
      m_a1 = RooRealVar("a1", "a1", Value, Minimum, Maximum);
    } else if(ParameterName == "a2") {
      m_a2 = RooRealVar("a2", "a2", Value, Minimum, Maximum);
    } else if(ParameterName == "b1") {
      m_b1 = RooRealVar("b1", "b1", Value, Minimum, Maximum);
    } else if(ParameterName == "b2") {
      m_b2 = RooRealVar("b2", "b2", Value, Minimum, Maximum);
    } else {
      continue;
    }
  }
}

void DeltaEFit::FitDeltaE(const std::string &Filename, const std::string &TagMode, bool DoBinnedFit, bool DoUnbinnedFit) {
  using namespace RooFit;
  TH1D h1("h1", "h1", 1000, 0.08, 0.08);
  m_Tree->Draw("DeltaE >> h1", "LuminosityWeight", "goff");
  RooDataHist BinnedData("BinnedData", "BinnedData", RooArgList(m_DeltaE), &h1);
  RooDataSet UnbinnedData("UnbinnedData", "UnbinnedData", m_Tree, RooArgList(m_DeltaE, m_LuminosityWeight), "", "LuminosityWeight");
  RooGaussian Gauss1("Gauss1", "Gauss1", m_DeltaE, m_Mean1, m_Sigma1);
  RooGaussian Gauss2("Gauss2", "Gauss2", m_DeltaE, m_Mean2, m_Sigma2);
  //RooPolynomial Polynomial("Polynomial", "Polynomial", m_DeltaE, RooArgList(m_a, m_b));
  RooGenericPdf Polynomial("Polynomial", "DeltaE < 0 ? 1 + a1*DeltaE + b1*DeltaE*DeltaE : 1 + a2*DeltaE + b2*DeltaE*DeltaE", RooArgList(m_DeltaE, m_a1, m_a2, m_b1, m_b2));
  RooAddPdf Model("Model", "Model", RooArgList(Gauss1, Gauss2, Polynomial), RooArgList(m_Nsig1, m_Nsig2, m_Nbkg));
  if(DoBinnedFit) {
    Model.fitTo(BinnedData);
    if(DoUnbinnedFit) {
      Model.fitTo(UnbinnedData);
    }
  }
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
  RooPlot *Frame = m_DeltaE.frame();
  Frame->SetTitle((TagMode + std::string(" Single Tag #Delta E fit;#Delta E (GeV);Events")).c_str());
  UnbinnedData.plotOn(Frame, Binning(100));
  Model.plotOn(Frame, LineColor(kBlue));
  RooHist *Pull = Frame->pullHist();
  Model.plotOn(Frame, Components(Polynomial), LineColor(kBlue), LineStyle(kDashed));
  Model.plotOn(Frame, Components(RooArgList(Gauss1, Gauss2)), LineColor(kRed));
  if(m_Tree->GetEntries("DataSetType == 0") == 0) {
    TH1D h2("h2", "h2", 100, -0.08, 0.08);
    m_Tree->Draw("DeltaE >> h2", "LuminosityWeight*(DataSetType == 1 && IsSameDMother == 1 && PIDTrue == 1)", "goff");
    RooDataHist TrueSignal("TrueSignal", "TrueSignal", RooArgList(m_DeltaE), &h2);
    TrueSignal.plotOn(Frame);
  }
  Frame->Draw();
  Pad2.cd();
  RooPlot *PullFrame = m_DeltaE.frame();
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

void DeltaEFit::SaveParameters(const std::string &Filename) const {
  std::ofstream OutputFile(Filename);
  OutputFile << "Nsig1 " << m_Nsig1.getValV() << "\n";
  OutputFile << "Nsig2 " << m_Nsig2.getValV() << "\n";
  OutputFile << "Nbkg " << m_Nbkg.getValV() << "\n";
  OutputFile << "Mean1 " << m_Mean1.getValV() << "\n";
  OutputFile << "Mean2 " << m_Mean2.getValV() << "\n";
  OutputFile << "Sigma1 " << m_Sigma1.getValV() << "\n";
  OutputFile << "Sigma2 " << m_Sigma2.getValV() << "\n";
  OutputFile << "a1 " << m_a1.getValV() << "\n";
  OutputFile << "b1 " << m_b1.getValV() << "\n";
  OutputFile << "a2 " << m_a1.getValV() << "\n";
  OutputFile << "b2 " << m_b1.getValV() << "\n";
  OutputFile.close();
}
