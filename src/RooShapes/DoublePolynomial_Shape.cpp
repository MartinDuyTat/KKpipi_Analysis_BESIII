// Martin Duy Tat 17th November 2021

#include<string>
#include<vector>
#include<stdexcept>
#include"RooRealVar.h"
#include"RooGenericPdf.h"
#include"RooArgList.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoublePolynomial_Shape.h"
#include"Utilities.h"
#include"Unique.h"

DoublePolynomial_Shape::DoublePolynomial_Shape(const std::string &Name, const Settings &settings, RooRealVar *x): FitShape(Name, settings, x), m_Order(m_Settings.getI("PolynomialOrder")) {
  if(m_Order <= 0) {
    throw std::invalid_argument("Polynomial order must be strictly positive");
  }
}

void DoublePolynomial_Shape::Initialize() {
  RooArgList PolynomialParameters(*m_x);
  // Left and right formulas for polynomials
  std::vector<std::string> PolynomialFormulas{"abs(1", "abs(1"};
  int Parameter_index = 1;
  for(int i = 1; i <= 2; i++) {
    for(int j = 0; j < m_Order; j++) {
      char CoefficientLabel = 'a' + j;
      PolynomialFormulas[i - 1] += " + ";
      m_Parameters.insert({CoefficientLabel + std::to_string(i), Utilities::load_param(m_Settings, m_Name + "_" + CoefficientLabel + std::to_string(i))});
      PolynomialParameters.add(*m_Parameters[CoefficientLabel + std::to_string(i)]);
      PolynomialFormulas[i - 1] += "@" + std::to_string(Parameter_index);
      Parameter_index++;
      for(int k = 0; k < j + 1; k++) {
	PolynomialFormulas[i - 1] += "*@0";
      }
      CoefficientLabel++;
    }
    PolynomialFormulas[i - 1] += ")";
  }
  std::string PolynomialFormula = "@0 < 0 ? " + PolynomialFormulas[0] + " : " + PolynomialFormulas[1];
  m_PDF = Unique::create<RooGenericPdf*>(m_Name + "_DoublePolynomial", PolynomialFormula.c_str(), PolynomialParameters);
}
