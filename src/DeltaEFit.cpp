// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include"TTree.h"
#include"TH1D.h"
#include"TCanvas.h"
#include"TPad.h"
#include"TLine.h"
#include"RooRealVar.h"
#include"RooDataHist.h"
#include"RooDataSet.h"
#include"RooAbsPdf.h"
#include"RooArgList.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"RooFitResult.h"
#include"DeltaEFit.h"
#include"DeltaEFitModel.h"
#include"Settings.h"

DeltaEFit::DeltaEFit(TTree *Tree, const Settings &settings):
                     m_Tree(Tree),
		     m_Settings(settings),
		     m_DeltaE(RooRealVar("DeltaE", "DeltaE", -0.08, 0.08)),
                     m_LuminosityWeight("LuminosityWeight", "LuminosityWeight", 1.0, 0.0, 10.0) {
  m_Tree->SetBranchStatus("*", 0);
  m_Tree->SetBranchStatus("DeltaE" , 1);
  m_Tree->SetBranchStatus("DataSetType", 1);
  m_Tree->SetBranchStatus("IsSameDMother", 1);
  m_Tree->SetBranchStatus("PIDTrue", 1);
  m_Tree->SetBranchStatus("LuminosityWeight", 1);
}

void DeltaEFit::FitDeltaE() {
  using namespace RooFit;
  DeltaEFitModel FitModel(m_Settings, &m_DeltaE);
  TH1D h1("h1", "h1", 1000, 0.08, 0.08);
  m_Tree->Draw("DeltaE >> h1", "LuminosityWeight", "goff");
  RooDataHist BinnedData("BinnedData", "BinnedData", RooArgList(m_DeltaE), &h1);
  RooDataSet UnbinnedData("UnbinnedData", "UnbinnedData", m_Tree, RooArgList(m_DeltaE, m_LuminosityWeight), "", "LuminosityWeight");
  RooFitResult *Results;
  auto Model = FitModel.GetModel();
  if(m_Settings.get("FitType") != "NoFit") {
    Results = Model->fitTo(BinnedData, Save());
    if(m_Settings.get("FitType") == "UnbinnedFit") {
      Results = Model->fitTo(UnbinnedData, Save());
    }
    SaveParameters(Results);
  }
  SavePlot(UnbinnedData, Model);
}

void DeltaEFit::SavePlot(const RooDataSet &UnbinnedData, RooAbsPdf *Model) const {
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
  RooPlot *Frame = m_DeltaE.frame();
  Frame->SetTitle((m_Settings.get("Mode") + std::string(" Single Tag #Delta E fit;#Delta E (GeV);Events")).c_str());
  UnbinnedData.plotOn(Frame, Binning(100));
  Model->plotOn(Frame, LineColor(kBlue));
  RooHist *Pull = Frame->pullHist();
  Model->plotOn(Frame, Components((m_Settings.get("Mode") + "_SingleTag_Combinatorial_" + m_Settings["Combinatorial"].get("CombinatorialShape")).c_str()), LineColor(kBlue), LineStyle(kDashed));
  Model->plotOn(Frame, Components((m_Settings.get("Mode") + "_SingleTag_Signal_DoubleGaussian").c_str()), LineColor(kRed));
  Frame->Draw();
  Pad2.cd();
  RooPlot *PullFrame = m_DeltaE.frame();
  PullFrame->addObject(Pull, "P E1");
  TLine *Line = new TLine(-0.08, 0.0, 0.08, 0.0);
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
  c1.SaveAs(m_Settings.get("PlotFilename").c_str());
}

void DeltaEFit::SaveParameters(RooFitResult *Results) const {
  Results->Print("V");
  std::ofstream OutputFile(m_Settings.get("ResultsFilename"));
  OutputFile << "status " << Results->status() << "\n";
  OutputFile << "covQual " << Results->covQual() << "\n";
  RooArgList floating_param = Results->floatParsFinal();
  for(int i = 0; i < floating_param.getSize(); i++) {
    RooRealVar *param = static_cast<RooRealVar*>(floating_param.at(i));
    OutputFile << param->GetName() << " " << param->getVal() << "\n";
    OutputFile << param->GetName() << "_err " << param->getError() << "\n";
  }
  OutputFile.close();
}
