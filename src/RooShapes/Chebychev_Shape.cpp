// Martin Duy Tat 18th November 2021

#include<stdexcept>
#include<string>
#include"RooRealVar.h"
#include"RooChebychev.h"
#include"RooArgList.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/Chebychev_Shape.h"
#include"Utilities.h"
#include"Unique.h"

Chebychev_Shape::Chebychev_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name + "_Combinatorial", settings, x), m_Order(m_Settings.getI("PolynomialOrder")) {
  if(m_Order <= 0) {
    throw std::invalid_argument("Polynomial order must be strictly positive");
  }
}

void Chebychev_Shape::Initialize() {
  RooArgList Coefficients;
  for(int i = 0; i < m_Order; i++) {
    m_Parameters.insert({"c" + std::to_string(i), Utilities::load_param(m_Settings, m_Name + "_c" + std::to_string(i))});
    Coefficients.add(*m_Parameters.at("c" + std::to_string(i)));
  }
  m_PDF = Unique::create<RooChebychev*>(m_Name + "_Chebychev", "", *m_x, Coefficients);
}
