// Martin Duy Tat 17th November 2021

#include"RooRealVar.h"
#include"RooGaussian.h"
#include"RooAddPdf.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"Utilities.h"
#include"Unique.h"

DoubleGaussian_Shape::DoubleGaussian_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x) {
}

void DoubleGaussian_Shape::Initialize() {
  m_Parameters.insert({"mu1", Utilities::load_param(m_Settings, m_Name + "_mu1")});
  m_Parameters.insert({"mu2", Utilities::load_param(m_Settings, m_Name + "_mu2")});
  m_Parameters.insert({"sigma1", Utilities::load_param(m_Settings, m_Name + "_sigma1")});
  m_Parameters.insert({"sigma2", Utilities::load_param(m_Settings, m_Name + "_sigma2")});
  m_Parameters.insert({"frac", Utilities::load_param(m_Settings, m_Name + "_frac")});
  auto Gaussian1 = Unique::create<RooGaussian*>(m_Name + "_Gaussian1", "", *m_x, *m_Parameters["mu1"], *m_Parameters["sigma1"]);
  auto Gaussian2 = Unique::create<RooGaussian*>(m_Name + "_Gaussian2", "", *m_x, *m_Parameters["mu2"], *m_Parameters["sigma2"]);
  m_PDF = Unique::create<RooAddPdf*>(m_Name + "_DoubleGaussian", "", *Gaussian1, *Gaussian2, *m_Parameters["frac"]);
}
