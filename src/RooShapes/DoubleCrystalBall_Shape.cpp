// Martin Duy Tat 17th November 2021

#include"RooRealVar.h"
#include"RooCBShape.h"
#include"RooAddPdf.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"
#include"Utilities.h"
#include"Unique.h"

DoubleCrystalBall_Shape::DoubleCrystalBall_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x) {
}

void DoubleCrystalBall_Shape::Initialize() {
  m_Parameters.insert({"mu1", Utilities::load_param(m_Settings, m_Name + "_mu1")});
  m_Parameters.insert({"mu2", Utilities::load_param(m_Settings, m_Name + "_mu2")});
  m_Parameters.insert({"sigma1", Utilities::load_param(m_Settings, m_Name + "_sigma1")});
  m_Parameters.insert({"sigma2", Utilities::load_param(m_Settings, m_Name + "_sigma2")});
  m_Parameters.insert({"alpha1", Utilities::load_param(m_Settings, m_Name + "_alpha1")});
  m_Parameters.insert({"alpha2", Utilities::load_param(m_Settings, m_Name + "_alpha2")});
  m_Parameters.insert({"n1", Utilities::load_param(m_Settings, m_Name + "_n1")});
  m_Parameters.insert({"n2", Utilities::load_param(m_Settings, m_Name + "_n2")});
  m_Parameters.insert({"frac", Utilities::load_param(m_Settings, m_Name + "_frac")});
  auto CrystalBall1 = Unique::create<RooCBShape*>(m_Name + "_CrystalBall1", "", *m_x, *m_Parameters["mu1"], *m_Parameters["sigma1"], *m_Parameters["alpha1"], *m_Parameters["n1"]);
  auto CrystalBall2 = Unique::create<RooCBShape*>(m_Name + "_CrystalBall2", "", *m_x, *m_Parameters["mu2"], *m_Parameters["sigma2"], *m_Parameters["alpha2"], *m_Parameters["n2"]);
  m_PDF = Unique::create<RooAddPdf*>(m_Name + "_DoubleCrystalBall", "", *CrystalBall1, *CrystalBall2, *m_Parameters["frac"]);
}
