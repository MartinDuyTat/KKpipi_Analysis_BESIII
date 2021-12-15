// Martin Duy Tat 26th November 2021

#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include"TTree.h"
#include"TPad.h"
#include"TCanvas.h"
#include"TAxis.h"
#include"TLine.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooSimultaneous.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"DoubleTagYield.h"
#include"Settings.h"
#include"BinnedDataLoader.h"
#include"BinnedFitModel.h"
#include"Category.h"

DoubleTagYield::DoubleTagYield(const Settings &settings, TTree *Tree): m_SignalMBC("SignalMBC", "", 1.83, 1.8865),
								       m_TagMBC("TagMBC", "", 1.83, 1.8865),
								       m_Settings(settings), m_Tree(Tree) {
  m_SignalMBC.setBins(5000, "cache");
}

void DoubleTagYield::DoFit() {
  using namespace RooFit;
  BinnedDataLoader DataLoader(m_Settings, m_Tree, &m_SignalMBC, &m_TagMBC);
  RooDataSet *DataSet = DataLoader.GetDataSet();
  BinnedFitModel FitModel(m_Settings, &m_SignalMBC, &m_TagMBC);
  RooSimultaneous *Model = FitModel.GetPDF();
  // Perform an initial fit
  auto Result = Model->fitTo(*DataSet, Save(), NumCPU(4));
  // Any bins with less than 0.5 combinatorial background events are set constant
  std::vector<std::string> Categories = DataLoader.GetCategoryObject()->GetCategories();
  for(const auto &Category : Categories) {
    RooRealVar *CombinatorialYield = static_cast<RooRealVar*>(FitModel.m_Yields[Category + "_CombinatorialYield"]);
    if(CombinatorialYield->getVal() < 0.5) {
      CombinatorialYield->setConstant();
    }
  }
  for(const auto &Parameter : FitModel.m_Parameters) {
    if(Parameter.first == "End") {
      continue;
    }
    Parameter.second->setConstant();
  }
  // Perform a second fit if fit is binned
  if(!(m_Settings.contains("Inclusive_fit") && m_Settings.getB("Inclusive_fit"))) {
    Result = Model->fitTo(*DataSet, Save(), NumCPU(4));
  }
  Result->Print();
  PlotProjections(&DataLoader, &FitModel);
  SaveSignalYields(FitModel, Result, *DataLoader.GetCategoryObject());
}

void DoubleTagYield::PlotProjections(BinnedDataLoader *DataLoader, BinnedFitModel *FitModel) {
  using namespace RooFit;
  auto Model = FitModel->GetPDF();
  Category *category = DataLoader->GetCategoryObject();
  RooCategory *CategoryVariable = category->GetCategoryVariable();
  RooDataSet *DataSet = DataLoader->GetDataSet();
  for(const auto &Category : category->GetCategories()) {
    int Bin = category->GetSignalBinNumber(Category);
    TCanvas c1((Category + "_c1").c_str(), "", 1600, 1200);
    TPad Pad1((Category + "_Pad1").c_str(), "", 0.0, 0.25, 1.0, 1.0);
    TPad Pad2((Category + "_Pad2").c_str(), "", 0.0, 0.0, 1.0, 0.25);
    Pad1.Draw();
    Pad2.Draw();
    Pad1.SetBottomMargin(0.1);
    Pad1.SetTopMargin(0.1);
    Pad1.SetBorderMode(0);
    Pad2.SetBorderMode(0);
    Pad2.SetBottomMargin(0.1);
    Pad2.SetTopMargin(0.05);
    Pad1.cd();
    RooPlot *Frame = m_SignalMBC.frame();
    std::string TagMode = m_Settings.get("Mode");
    Frame->SetTitle((TagMode + " Double Tag M_{BC} Bin " + std::to_string(Bin) + "; M_{BC} (GeV); Events").c_str());
    DataSet->plotOn(Frame, Binning(m_Settings.getI("Bins_in_plots")), Cut((std::string(CategoryVariable->GetName()) + "==" + std::string(CategoryVariable->GetName()) + "::" + Category).c_str()));
    Model->plotOn(Frame, LineColor(kBlue), Slice(*CategoryVariable, Category.c_str()), ProjWData(*CategoryVariable, *DataSet));
    RooHist *Pull = Frame->pullHist();
    Model->plotOn(Frame, LineColor(kBlue), Components("Argus"), LineStyle(kDashed), Slice(*CategoryVariable, Category.c_str()), ProjWData(*CategoryVariable, *DataSet));
    if(FitModel->m_PeakingBackgroundShapes.size() > 0) {
      std::string PeakingList = "";
      for(auto iter = FitModel->m_PeakingBackgroundShapes.begin(); iter != FitModel->m_PeakingBackgroundShapes.end(); iter++) {
	if(iter != FitModel->m_PeakingBackgroundShapes.begin()) {
	  PeakingList += std::string(",");
	}
	PeakingList += iter->second->GetPDF()->GetName();
      }
      Model->plotOn(Frame, LineColor(kMagenta), Slice(*CategoryVariable, Category.c_str()), ProjWData(*CategoryVariable, *DataSet), Components(PeakingList.c_str()));
    }
    Frame->Draw();
    Pad2.cd();
    RooPlot *PullFrame = m_SignalMBC.frame();
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
    std::string PlotFilename = m_Settings.get("MBCPlotFilenamePrefix") + "_" + Category + ".png";
    c1.SaveAs(PlotFilename.c_str());
  }
}

void DoubleTagYield::SaveSignalYields(const BinnedFitModel &FitModel, RooFitResult *Result, const Category &category) const {
  double Fraction = FitModel.GetFractionInSignalRegion();
  std::ofstream Outfile(m_Settings.get("FittedSignalYieldsFile"));
  Outfile << std::fixed << std::setprecision(4);
  Outfile << "* KKpipi vs " << m_Settings.get("Mode") << " double tag yield fit results\n\n";
  Outfile << "status " << Result->status() << "\n";
  Outfile << "covQual " << Result->covQual() << "\n\n";
  Outfile << "* Fraction of events inside signal window: " << Fraction << "\n\n";
  for(const auto & cat : category.GetCategories()) {
    std::string Name = cat + "_SignalYield";
    Outfile << Name << "     " << std::setw(8) << std::right << FitModel.m_Yields.at(Name)->getVal()*Fraction << "\n";
    Outfile << Name << "_err " << std::setw(8) << std::right << static_cast<RooRealVar*>(FitModel.m_Yields.at(Name))->getError()*Fraction << "\n";
  }
  Outfile.close();
}
