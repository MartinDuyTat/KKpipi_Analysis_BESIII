// Martin Duy Tat 17th November 2021

#include"RooRealVar.h"
#include"RooGenericPdf.h"
#include"RooArgList.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoublePolynomial_Shape.h"
#include"Utilities.h"
#include"Unique.h"

DoublePolynomial_Shape::DoublePolynomial_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x) {
}

void DoublePolynomial_Shape::Initialize() {
  m_Parameters.insert({"a1", Utilities::load_param(m_Settings, m_Name + "_a1")});
  m_Parameters.insert({"a2", Utilities::load_param(m_Settings, m_Name + "_a2")});
  m_Parameters.insert({"b1", Utilities::load_param(m_Settings, m_Name + "_b1")});
  m_Parameters.insert({"b2", Utilities::load_param(m_Settings, m_Name + "_b2")});
  m_PDF = Unique::create<RooGenericPdf*>(m_Name + "_DoublePolynomial", "@0 < 0 ? abs(1 + @1*@0 + @2*@0*@0) : abs(1 + @3*@0 + @4*@0*@0)", RooArgList(*m_x, *m_Parameters["a1"], *m_Parameters["b1"], *m_Parameters["a2"], *m_Parameters["b2"]));
}
