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
#include"TCut.h"
#include"TFile.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooSimultaneous.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"RooMsgService.h"
#include"RooStats/SPlot.h"
#include"RooStats/RooStatsUtils.h"
#include"DoubleTagYield.h"
#include"Settings.h"
#include"BinnedDataLoader.h"
#include"BinnedFitModel.h"
#include"Category.h"

DoubleTagYield::DoubleTagYield(const Settings &settings, TTree *Tree): m_SignalMBC("SignalMBC", "", 1.83, 1.8865),
								       m_Settings(settings), m_Tree(Tree) {
  for(int i = 0; i < 2; i++) {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Plotting);
  }
  if(!m_Settings.getB("FullyReconstructed")) {
    m_SignalMBC = RooRealVar(m_Settings.get("FitVariable").c_str(), "", m_Settings.getD("FitRange_low"), m_Settings.getD("FitRange_high"));
  }
  m_SignalMBC.setBins(200, "cache");
}

void DoubleTagYield::DoFit() {
  using namespace RooFit;
  BinnedDataLoader DataLoader(m_Settings, m_Tree, &m_SignalMBC);
  RooDataSet *DataSet = DataLoader.GetDataSet();
  BinnedFitModel FitModel(m_Settings, &m_SignalMBC);
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
  if(Categories.size() > 1) {
    Result = Model->fitTo(*DataSet, Save(), NumCPU(4));
  }
  Result->Print();
  PlotProjections(&DataLoader, &FitModel);
  SaveSignalYields(FitModel, Result, *DataLoader.GetCategoryObject());
  if(m_Settings.contains("sPlotReweight") && m_Settings.getB("sPlotReweight")) {
    sPlotReweight(*DataSet, FitModel);
  }
}

void DoubleTagYield::PlotProjections(BinnedDataLoader *DataLoader, BinnedFitModel *FitModel) {
  using namespace RooFit;
  auto Model = FitModel->GetPDF();
  Category *category = DataLoader->GetCategoryObject();
  RooCategory *CategoryVariable = category->GetCategoryVariable();
  RooDataSet *DataSet = DataLoader->GetDataSet();
  for(const auto &Category : category->GetCategories()) {
    int SignalBin = category->GetSignalBinNumber(Category);
    int TagBin = category->GetTagBinNumber(Category);
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
    std::string Title = TagMode + " Double Tag M_{BC}";
    if(SignalBin != 0) {
      Title += ", KK#pi#pi bin " + std::to_string(SignalBin);
    } else {
      Title += ", inclusive KK#pi#pi phase space";
    }
    if(TagMode == "KSpipi" || TagMode == "KSKK" || TagMode == "KKpipi") {
      Title += ", tag bin " + std::to_string(TagBin);
    }
    if(TagMode.substr(0, 2) == "KL") {
      Title += "; M_{miss}^{2} (GeV^{2}); Events";
    } else if(TagMode == "KeNu") {
      Title += "; U_{miss} (GeV); Events";
    } else {
      Title += "; M_{BC} (GeV); Events";
    }
    Frame->SetTitle(Title.c_str());
    DataSet->plotOn(Frame, Binning(m_Settings.getI("Bins_in_plots")), Cut((std::string(CategoryVariable->GetName()) + "==" + std::string(CategoryVariable->GetName()) + "::" + Category).c_str()));
    Model->plotOn(Frame, LineColor(kBlue), Slice(*CategoryVariable, Category.c_str()), ProjWData(*CategoryVariable, *DataSet));
    RooHist *Pull = Frame->pullHist();
    Model->plotOn(Frame, LineColor(kBlue), Components("Combinatorial*"), LineStyle(kDashed), Slice(*CategoryVariable, Category.c_str()), ProjWData(*CategoryVariable, *DataSet));
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
    TLine *Line = new TLine(m_SignalMBC.getMin(), 0.0, m_SignalMBC.getMax(), 0.0);
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
  std::ofstream Outfile(m_Settings.get("FittedSignalYieldsFile"));
  Outfile << std::fixed << std::setprecision(4);
  Outfile << "* KKpipi vs " << m_Settings.get("Mode") << " double tag yield fit results\n\n";
  Outfile << "status " << Result->status() << "\n";
  Outfile << "covQual " << Result->covQual() << "\n\n";
  if(m_Settings.getB("FullyReconstructed")) {
    Outfile << "* These yields are after the sideband has been subtracted off the signal yield\n\n";
  }
  for(const auto & cat : category.GetCategories()) {
    std::string Name = cat + "_SignalYield";
    double Sideband = 0.0;
    if(m_Settings.getB("FullyReconstructed")) {
      Sideband += GetSidebandYield(category.GetSignalBinNumber(cat), category.GetTagBinNumber(cat));
    }
    Outfile << Name << "          " << std::setw(8) << std::right << FitModel.m_Yields.at(Name)->getVal() - Sideband << "\n";
    Outfile << Name << "_err      " << std::setw(8) << std::right << static_cast<RooRealVar*>(FitModel.m_Yields.at(Name))->getError() << "\n";
    if(m_Settings.getB("FullyReconstructed")) {
      Outfile << Name << "_sideband " << std::setw(8) << std::right << Sideband << "\n";
    }
  }
  Outfile.close();
}

double DoubleTagYield::GetSidebandYield(int SignalBin, int TagBin) const {
  TCut SidebandCut("TagMBC > 1.84 && TagMBC < 1.85 && SignalMBC > 1.86 && SignalMBC < 1.87");
  if(SignalBin != 0) {
    SidebandCut = SidebandCut && TCut(("SignalBin == " + std::to_string(SignalBin)).c_str());
  }
  if(TagBin != 0) {
    SidebandCut = SidebandCut && TCut(("TagBin == " + std::to_string(TagBin)).c_str());
  }
  return static_cast<double>(m_Tree->GetEntries(SidebandCut));
}

void DoubleTagYield::sPlotReweight(RooDataSet &Data, BinnedFitModel &FitModel) {
  auto FloatingParameters = FitModel.GetPDF()->getParameters(Data);
  for(auto Parameter : *FloatingParameters) {
    auto CastedParameter = dynamic_cast<RooRealVar*>(Parameter);
    if(CastedParameter) {
      CastedParameter->setConstant(true);
    }
  }
  RooArgList YieldParameters;
  for(const auto &CategoryString : FitModel.m_Category.GetCategories()) {
    auto YieldParameter = static_cast<RooRealVar*>(FitModel.m_Yields.at(CategoryString + "_SignalYield"));
    YieldParameter->setConstant(false);
    YieldParameters.add(*YieldParameter);
  }
  RooStats::SPlot("sData", "", Data, FitModel.GetPDF(), YieldParameters);
  TFile Outfile(m_Settings.get("sPlotFilename").c_str(), "RECREATE");
  auto Tree = RooStats::GetAsTTree(m_Settings.get("TreeName").c_str(), m_Settings.get("TreeName").c_str(), Data);
  Tree->Write();
  Outfile.Close();
}
