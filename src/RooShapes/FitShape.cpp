// Martin Duy Tat 17th November 2021

#include<string>
#include"RooAbsPdf.h"
#include"RooRealVar.h"
#include"RooShapes/FitShape.h"
#include"Settings.h"
#include"Utilities.h"

FitShape::FitShape(const std::string &Name, const Settings &settings, RooRealVar *x):  m_Settings(settings), m_Name(Name), m_x(x) {
  m_Yield = Utilities::load_param(m_Settings, m_Name + "_yield");
}

FitShape::~FitShape() {
}

RooAbsPdf* FitShape::GetPDF() {
  if(!m_PDF) {
    Initialize();
  }
  return m_PDF;
}

RooRealVar* FitShape::GetYield() {
  if(!m_Yield) {
    Initialize();
  }
  return m_Yield;
}
