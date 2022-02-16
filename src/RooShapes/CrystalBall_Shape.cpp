// Martin Duy Tat 17th November 2021

#include"RooRealVar.h"
#include"RooCBShape.h"
#include"RooAddPdf.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/CrystalBall_Shape.h"
#include"Utilities.h"
#include"Unique.h"

CrystalBall_Shape::CrystalBall_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x) {
}

void CrystalBall_Shape::Initialize() {
  m_Parameters.insert({"mu", Utilities::load_param(m_Settings, m_Name + "_mu")});
  m_Parameters.insert({"sigma", Utilities::load_param(m_Settings, m_Name + "_sigma")});
  m_Parameters.insert({"alpha", Utilities::load_param(m_Settings, m_Name + "_alpha")});
  m_Parameters.insert({"n", Utilities::load_param(m_Settings, m_Name + "_n")});
  m_PDF = Unique::create<RooCBShape*>(m_Name + "_CrystalBall", "", *m_x, *m_Parameters["mu"], *m_Parameters["sigma"], *m_Parameters["alpha"], *m_Parameters["n"]);
}
