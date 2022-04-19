// Martin Duy Tat 17th November 2021

#include<string>
#include"RooAbsPdf.h"
#include"RooRealVar.h"
#include"RooFormulaVar.h"
#include"RooArgSet.h"
#include"Settings.h"
#include"Utilities.h"
#include"Unique.h"
#include"RooShapes/FitShape.h"

FitShape::FitShape(const std::string &Name, const Settings &settings, RooRealVar *x):  m_Settings(settings), m_Name(Name), m_x(x) {
  if(m_Settings.contains(m_Name + "_yield")) {
    m_Yield = Utilities::load_param(m_Settings, m_Name + "_yield");
  } else {
    m_Yield = Unique::create<RooRealVar*>((m_Name + "_yield").c_str(), "", 0.0, 0.0, 1000000.0);
  }
}

FitShape::~FitShape() {
  if(m_Yield) {
    delete m_Yield;
  }
}

RooAbsPdf* FitShape::GetPDF() {
  if(!m_PDF) {
    Initialize();
  }
  return m_PDF;
}

RooAbsReal* FitShape::GetYield() {
  if(!m_Yield) {
    Initialize();
  }
  return m_Yield;
}

void FitShape::UseRelativeYield(RooRealVar *SignalYield, double BackgroundToSignalYieldRatio) {
  if(!m_Yield) {
    delete m_Yield;
  }
  auto BkgToSigRatioVar = Unique::create<RooRealVar*>(m_Name + "_BackgroundToSignalYieldRatio", "", BackgroundToSignalYieldRatio);
  m_Yield = Unique::create<RooFormulaVar*>(m_Name + "_yield", "@0*@1", RooArgSet(*SignalYield, *BkgToSigRatioVar));
}
