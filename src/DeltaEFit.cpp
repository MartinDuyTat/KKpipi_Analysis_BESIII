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
#include"RooPolynomial.h"
#include"RooAddPdf.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"DeltaEFit.h"

DeltaEFit::DeltaEFit(TTree *Tree):
                     m_Tree(Tree),
		     m_DeltaE(RooRealVar("DeltaE", "DeltaE", -0.08, 0.08)),
                     m_Nsig1(RooRealVar("Nsig1", "Nsig1", m_Tree->GetEntries()/2, 0, m_Tree->GetEntries())),
                     m_Nsig2(RooRealVar("Nsig2", "Nsig2", m_Tree->GetEntries()/3, 0, m_Tree->GetEntries())),
		     m_Nbkg(RooRealVar("Nbkg", "Nbkg", m_Tree->GetEntries()/2, 0, m_Tree->GetEntries())),
		     m_Mean1(RooRealVar("Mean1", "Mean1", 0.0, -0.02, 0.02)),
		     m_Mean2(RooRealVar("Mean2", "Mean2", 0.0, -0.02, 0.02)),
		     m_Sigma1(RooRealVar("Sigma1", "Sigma1", 0.00001, 0.03)),
                     m_Sigma2(RooRealVar("Sigma2", "Sigma2", 0.00001, 0.03)),
                     m_a(RooRealVar("a", "a", 0.0, -1000.0, 1000.0)),
                     m_b(RooRealVar("b", "b", 0.0, -1000.0, 1000.0)) {
  m_Tree->SetBranchStatus("*", 0);
  m_Tree->SetBranchStatus("DeltaE" , 1);
}

void DeltaEFit::FitDeltaE(const std::string &Filename, const std::string &TagMode, bool DoUnbinnedFit) {
  using namespace RooFit;
  TH1D h1("h1", "h1", 1000, 0.08, 0.08);
  m_Tree->Draw("DeltaE >> h1", "", "goff");
  m_Sigma1.setVal(h1.GetStdDev());
  m_Sigma2.setVal(h1.GetStdDev());
  RooDataHist BinnedData("BinnedData", "BinnedData", RooArgList(m_DeltaE), &h1);
  RooDataSet UnbinnedData("UnbinnedData", "UnbinnedData", m_Tree, RooArgList(m_DeltaE));
  RooGaussian Gauss1("Gauss1", "Gauss1", m_DeltaE, m_Mean1, m_Sigma1);
  RooGaussian Gauss2("Gauss2", "Gauss2", m_DeltaE, m_Mean2, m_Sigma2);
  RooPolynomial Polynomial("Polynomial", "Polynomial", m_DeltaE, RooArgList(m_a, m_b));
  RooAddPdf Model("Model", "Model", RooArgList(Gauss1, Gauss2, Polynomial), RooArgList(m_Nsig1, m_Nsig2, m_Nbkg));
  Model.fitTo(BinnedData);
  if(DoUnbinnedFit) {
    Model.fitTo(UnbinnedData);
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
  OutputFile << "a " << m_a.getValV() << "\n";
  OutputFile << "b " << m_b.getValV() << "\n";
  OutputFile.close();
}
