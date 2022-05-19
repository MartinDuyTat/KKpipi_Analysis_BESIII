// Martin Duy Tat 17th November 2021

#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooAddPdf.h"
#include"RooArgList.h"
#include"RooFormulaVar.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussianRatio_Shape.h"
#include"Utilities.h"
#include"Unique.h"

DoubleGaussianRatio_Shape::DoubleGaussianRatio_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x) {
}

void DoubleGaussianRatio_Shape::Initialize() {
  m_Parameters.insert({"mu", Utilities::load_param(m_Settings, m_Name + "_mu")});
  m_Parameters.insert({"mu_f", Utilities::load_param(m_Settings, m_Name + "_mu_f")});
  m_Parameters.insert({"sigma", Utilities::load_param(m_Settings, m_Name + "_sigma")});
  m_Parameters.insert({"sigma_f", Utilities::load_param(m_Settings, m_Name + "_sigma_f")});
  m_Parameters.insert({"frac", Utilities::load_param(m_Settings, m_Name + "_frac")});
  auto mu2 = Unique::create<RooFormulaVar*>(m_Name + "_mu2", "@0*@1", RooArgList(*m_Parameters["mu"], *m_Parameters["mu_f"]));
  auto sigma2 = Unique::create<RooFormulaVar*>(m_Name + "_sigma2", "@0*@1", RooArgList(*m_Parameters["sigma"], *m_Parameters["sigma_f"]));
  auto Gaussian1 = Unique::create<RooGaussian*>(m_Name + "_Gaussian1", "", *m_x, *m_Parameters["mu"], *m_Parameters["sigma"]);
  auto Gaussian2 = Unique::create<RooGaussian*>(m_Name + "_Gaussian2", "", *m_x, *mu2, *sigma2);
  m_PDF = Unique::create<RooAddPdf*>(m_Name + "_DoubleGaussian", "", *Gaussian1, *Gaussian2, *m_Parameters["frac"]);
}
