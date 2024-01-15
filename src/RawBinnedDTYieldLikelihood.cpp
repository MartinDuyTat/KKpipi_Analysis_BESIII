// Martin Duy Tat 12th May 2023

#include<vector>
#include<string>
#include"TFile.h"
#include"RooAbsReal.h"
#include"RooArgSet.h"
#include"RooMsgService.h"
#include"RooAbsPdf.h"
#include"RooWorkspace.h"
#include"RawBinnedDTYieldLikelihood.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedDTYieldLikelihood::RawBinnedDTYieldLikelihood(const std::string &Tag,
						       const Settings &settings,
						       const std::string &TagCategory,
						       int ToyNumber):
  m_TagMode(Tag),
  m_TagCategory(TagCategory),
  m_Workspace(GetWorkspace(Tag, settings, ToyNumber)),
  m_FullLikelihood(std::move(GetFullLikelihood())),
  m_Variables(m_FullLikelihood->getVariables()),
  m_Order(GetYieldOrder(settings["BinningScheme"].getI("NumberBins"))) {
}

RawBinnedDTYieldLikelihood::RawBinnedDTYieldLikelihood(const std::string &Tag,
						       const Settings &settings,
						       const std::string &TagCategory,
						       const std::string &ToyName,
						       int ToyNumber):
  m_TagMode(Tag),
  m_TagCategory(TagCategory),
  m_Order(GetYieldOrder(settings["BinningScheme"].getI("NumberBins"))) {
  for(int i = 0; i < 2; i++) {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Fitting);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Minimization);
  }
  std::string FCFilename = settings.get("FeldmanCousinsToyPath");
  FCFilename = Utilities::ReplaceString(FCFilename, "TAG", Tag);
  FCFilename += "/" + ToyName + std::to_string(ToyNumber) + ".root";
  TFile File(FCFilename.c_str(), "READ");
  m_Workspace.reset(File.Get<RooWorkspace>("Workspace"));
  m_FullLikelihood = std::move(GetFullLikelihood());
  File.Close();
}

double RawBinnedDTYieldLikelihood::GetLogLikelihood(
  const std::vector<double> &PredictedBinYields) const {
  for(std::size_t i = 0; i < PredictedBinYields.size(); i++) {
    auto YieldVar = m_Variables->find(m_Order[i].c_str());
    static_cast<RooRealVar*>(YieldVar)->setVal(PredictedBinYields[i]);
  }
  double LogLikelihood = m_FullLikelihood->getVal();
  return LogLikelihood;
}

std::unique_ptr<RooWorkspace> RawBinnedDTYieldLikelihood::GetWorkspace(
  const std::string &Tag,
  const Settings &settings,
  int ToyNumber) const {
  auto DTYieldFilename = RawBinnedDTYields::GetFilename(Tag, settings, ToyNumber);
  DTYieldFilename = Utilities::ReplaceString(DTYieldFilename, ".txt", ".root");
  TFile File(DTYieldFilename.c_str(), "READ");
  std::unique_ptr<RooWorkspace> w{File.Get<RooWorkspace>("Workspace")};
  return w;
}

std::unique_ptr<RooAbsReal> RawBinnedDTYieldLikelihood::GetFullLikelihood() const {
  for(int i = 0; i < 2; i++) {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Fitting);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Minimization);
  }
  const std::string PDFName = "Simultaneous_PDF_KKpipi_vs_" + m_TagMode;
  auto Model = m_Workspace->pdf(PDFName.c_str());
  auto Data = m_Workspace->data("hmaster");
  if(!Data) {
    Data = m_Workspace->data("InputData");
  }
  return std::unique_ptr<RooAbsReal>{Model->createNLL(*Data)};
}

std::vector<std::string> RawBinnedDTYieldLikelihood::GetYieldOrder(
  int NumberBins) const {
  std::vector<std::string> Order;
  if(m_TagCategory == "CP") {
    std::string Prefix = "DoubleTag_CP_KKpipi_vs_" + m_TagMode + "_SignalBin";
    for(int Bin = 1; Bin <= NumberBins; Bin++) {
      std::string Label = Prefix + std::to_string(Bin) + "_SignalYield";
      Order.push_back(Label);
    }
  } else if(m_TagCategory == "Flavour") {
    std::string Prefix = "DoubleTag_Flavour_KKpipi_vs_";
    Prefix += m_TagMode + "_SignalBin";
    for(int Bin = -NumberBins; Bin <= NumberBins; Bin++) {
      if(Bin == 0) {
	continue;
      }
      std::string Label = Prefix + (Bin > 0 ? "P" : "M");
      Label += std::to_string(TMath::Abs(Bin)) + "_SignalYield";
      Order.push_back(Label);
    }
  } else if(m_TagCategory == "SCMB") {
    std::string Prefix = "DoubleTag_SCMB_KKpipi_vs_";
    Prefix += m_TagMode + "_SignalBin";
    for(int TagBin = 1; TagBin <= 8; TagBin++) {
      for(int SignalBin = -NumberBins; SignalBin <= NumberBins; SignalBin++) {
	if(SignalBin == 0) {
	  continue;
	}
	std::string Label = Prefix + (SignalBin > 0 ? "P" : "M");
	Label += std::to_string(TMath::Abs(SignalBin));
	Label += "_TagBin" + std::to_string(TagBin) + "_SignalYield";
	Order.push_back(Label);
      }
    }
  }
  return Order;
}
