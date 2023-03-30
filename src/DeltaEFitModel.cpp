// Martin Duy Tat 17th November 2021

#include<string>
#include<stdexcept>
#include<memory>
#include"RooAbsPdf.h"
#include"RooRealVar.h"
#include"DeltaEFitModel.h"
#include"Settings.h"
#include"Unique.h"
#include"Utilities.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussianRatio_Shape.h"
#include"RooShapes/DoublePolynomial_Shape.h"
#include"RooShapes/Chebychev_Shape.h"

DeltaEFitModel::DeltaEFitModel(const Settings &settings, RooRealVar *x): m_FullModel(nullptr), m_Settings(settings), m_x(x) {
  InitializeSignal();
  InitializeCombinatorial();
  InitializeFullShape();
}

DeltaEFitModel::~DeltaEFitModel() {
  delete m_FullModel;
}

RooAbsPdf* DeltaEFitModel::GetModel() const {
  return m_FullModel;
}

RooAbsPdf* DeltaEFitModel::GetModelComponent(int i) const {
  return m_ModelComponents[i]->GetPDF();
}

void DeltaEFitModel::InitializeSignal() {
  std::string Name = m_Settings.get("Mode") + "_SingleTag_Signal";
  m_ModelComponents.push_back(std::make_unique<DoubleGaussianRatio_Shape>(Name, m_Settings["Signal"], m_x));
  m_ModelPDFs.add(*m_ModelComponents.back()->GetPDF());
}

void DeltaEFitModel::InitializeCombinatorial() {
  std::string Name = m_Settings.get("Mode") + "_SingleTag_Combinatorial";
  std::string CombinatorialShape = m_Settings["Combinatorial"].get("CombinatorialShape");
  if(CombinatorialShape == "DoublePolynomial") {
    m_ModelComponents.push_back(std::make_unique<DoublePolynomial_Shape>(Name, m_Settings["Combinatorial"], m_x));
  } else if(CombinatorialShape == "Chebychev") {
    m_ModelComponents.push_back(std::make_unique<Chebychev_Shape>(Name, m_Settings["Combinatorial"], m_x));
  } else {
    throw std::invalid_argument("Unknown combinatorial shape");
  }
  m_ModelPDFs.add(*m_ModelComponents.back()->GetPDF());
}

void DeltaEFitModel::InitializeFullShape() {
  auto YieldRatio = Utilities::load_param(m_Settings["Signal"], m_Settings.get("Mode") + "_SingleTag_Signal_YieldFrac");
  m_FullModel = Unique::create<RooAddPdf*>((m_Settings.get("Mode") + "_SingleTag_Model").c_str(), "", m_ModelPDFs, *YieldRatio);
}
