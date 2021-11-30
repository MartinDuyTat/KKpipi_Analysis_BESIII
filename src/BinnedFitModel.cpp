// Martin Duy Tat 25th November 2021

#include<string>
#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooKeysPdf.h"
#include"RooArgusBG.h"
#include"RooAddPdf.h"
#include"RooSimultaneous.h"
#include"RooDataSet.h"
#include"TTree.h"
#include"BinnedFitModel.h"
#include"Settings.h"
#include"Utilities.h"
#include"Unique.h"

BinnedFitModel::BinnedFitModel(const Settings &settings, TTree *Tree, RooRealVar *SignalMBC): m_SignalMBC(SignalMBC),
											      m_Category(settings),
											      m_Settings(settings) {
  // First initialize the simultaneous fit with the correct categories for all bins
  auto CategoryVariable = m_Category.GetCategoryVariable();
  m_PDF = new RooSimultaneous(("Simultaneous_PDF_KKpipi_vs_" + settings.get("Mode")).c_str(), "", *CategoryVariable);
  // Set up the signal shape using signal MC
  InitializeSignalShape(Tree);
  // Then set up the combinatorial shape using an Argus function
  InitializeArgusShape();
  // Finally set up the simultaneous PDF for all bins
  InitializePDF();
}

BinnedFitModel::~BinnedFitModel() {
  if(m_PDF) {
    delete m_PDF;
  }
}

RooSimultaneous* BinnedFitModel::GetPDF() {
  return m_PDF;
}

void BinnedFitModel::InitializeSignalShape(TTree *Tree) {
  //auto frac = Unique::create<RooRealVar*>("frac", "", 0.4, 0.0, 1.0);
  m_Parameters.insert({"Mean1", Unique::create<RooRealVar*>("Mean1", "", 0.0, -0.001, 0.001)});
  //auto Mean2 = Unique::create<RooRealVar*>("Mean2", "", 0.0, -0.001, 0.001);
  m_Parameters.insert({"Sigma1", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_Sigma1")});
  //auto Sigma2 = Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_Sigma2");
  auto /*Gaussian1*/ Resolution = Unique::create<RooGaussian*>("Gaussian1", "", *m_SignalMBC, *m_Parameters["Mean1"], *m_Parameters["Sigma1"]);
  //auto Gaussian2 = Unique::create<RooGaussian*>("Gaussian2", "", *m_SignalMBC, *Mean2, *Sigma2);
  //auto Resolution = Unique::create<RooAddPdf*>("Resolution", "", RooArgList(*Gaussian1, *Gaussian2), *frac);
  RooDataSet MCSignal("MCSignal", "", Tree, RooArgList(*m_SignalMBC));
  auto SignalShape = Unique::create<RooKeysPdf*>("SignalShape", "", *m_SignalMBC, MCSignal);
  m_SignalShapeConv = Unique::create<RooFFTConvPdf*>("SignalShapeConv", "", *m_SignalMBC, *SignalShape, *Resolution);
}

void BinnedFitModel::InitializeArgusShape() {
  m_Parameters.insert({"End", Unique::create<RooRealVar*>("End", "", 1.8865)});
  m_Parameters.insert({"c", Utilities::load_param(m_Settings["MBC_Shape"], m_Settings.get("Mode") + "_DoubleTag_c")});
  m_Argus = Unique::create<RooArgusBG*>("Argus", "", *m_SignalMBC, *m_Parameters["End"], *m_Parameters["c"]);
}

RooAddPdf* BinnedFitModel::CreateBinPDF(const std::string &CategoryString) {
  m_Yields.insert({CategoryString + "_CombinatorialYield", Unique::create<RooRealVar*>((CategoryString + "_CombinatorialYield").c_str(), "", 1.0, 0.0, 1000.0)});
  m_Yields.insert({CategoryString + "_SignalYield", Unique::create<RooRealVar*>((CategoryString + "_SignalYield").c_str(), "", 10.0, 0.0, 1000.0)});
  return Unique::create<RooAddPdf*>((CategoryString + "_PDF").c_str(), "", RooArgList(*m_SignalShapeConv, *m_Argus), RooArgList(*m_Yields.at(CategoryString + "_SignalYield"), *m_Yields.at(CategoryString + "_CombinatorialYield")));
}

void BinnedFitModel::InitializePDF() {
  for(const auto &CategoryString : m_Category.GetCategories()) {
    RooAddPdf* BinPDF = CreateBinPDF(CategoryString);
    m_PDF->addPdf(*BinPDF, CategoryString.c_str());
  }
}

double BinnedFitModel::GetFractionInSignalRegion() const {
  using namespace RooFit;
  m_SignalMBC->setRange("SignalRange", 1.86, 1.87);
  return m_SignalShapeConv->createIntegral(*m_SignalMBC, NormSet(*m_SignalMBC), Range("SignalRange"))->getVal();
}
