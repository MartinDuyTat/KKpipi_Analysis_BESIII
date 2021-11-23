// Martin Duy Tat 17th November 2021

#include<string>
#include<stdexcept>
#include"RooAbsPdf.h"
#include"RooRealVar.h"
#include"DeltaEFitModel.h"
#include"Settings.h"
#include"Unique.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoublePolynomial_Shape.h"
#include"RooShapes/Chebychev_Shape.h"

DeltaEFitModel::DeltaEFitModel(const Settings &settings, RooRealVar *x): m_Settings(settings), m_x(x) {
  InitializeSignal();
  InitializeCombinatorial();
  InitializeFullShape();
}

DeltaEFitModel::~DeltaEFitModel() {
  for(auto Component : m_ModelComponents) {
    delete Component;
  }
  m_ModelComponents.clear();
}

RooAbsPdf* DeltaEFitModel::GetModel() const {
  return m_FullModel;
}

void DeltaEFitModel::InitializeSignal() {
  std::string Name = m_Settings.get("Mode") + "_SingleTag_Signal";
  m_ModelComponents.push_back(new DoubleGaussian_Shape(Name, m_Settings["Signal"], m_x));
  m_ModelPDFs.add(*m_ModelComponents.back()->GetPDF());
  m_ModelYields.add(*m_ModelComponents.back()->GetYield());
}

void DeltaEFitModel::InitializeCombinatorial() {
  std::string Name = m_Settings.get("Mode") + "_SingleTag_Combinatorial";
  std::string CombinatorialShape = m_Settings["Combinatorial"].get("CombinatorialShape");
  if(CombinatorialShape == "DoublePolynomial") {
    m_ModelComponents.push_back(new DoublePolynomial_Shape(Name, m_Settings["Combinatorial"], m_x));
  } else if(CombinatorialShape == "Chebychev") {
    m_ModelComponents.push_back(new Chebychev_Shape(Name, m_Settings["Combinatorial"], m_x));
  } else {
    throw std::invalid_argument("Unknown combinatorial shape");
  }
  m_ModelPDFs.add(*m_ModelComponents.back()->GetPDF());
  m_ModelYields.add(*m_ModelComponents.back()->GetYield());
}

void DeltaEFitModel::InitializeFullShape() {
  m_FullModel = Unique::create<RooAddPdf*>((m_Settings.get("Mode") + "_SingleTag_Model").c_str(), "", m_ModelPDFs, m_ModelYields);
}
