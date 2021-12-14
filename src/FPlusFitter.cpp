// Martin Duy Tat 14th December 2021

#include"TString.h"
#include"TMatrixTSym.h"
#include"RooRealVar.h"
#include"RooFormulaVar.h"
#include"RooArgList.h"
#include"RooMultiVarGaussian.h"
#include"RooFitResult.h"
#include"RooDataSet.h"
#include"TMath.h"
#include"Settings.h"
#include"Unique.h"
#include"FPlusFitter.h"

FPlusFitter::FPlusFitter(const Settings &settings): m_FPlus("FPlus", "", 0.0, 1.0), m_Settings(settings) {
}

void FPlusFitter::AddTag(const std::string &TagMode) {
  AddMeasurement(TagMode);
  AddPrediction(TagMode);
}

void FPlusFitter::InitializeAndFit() {
  // Set up trivial (diagonal, assuming all tags are independent) covariance matrix
  auto NumberTags = m_Uncertainties.size();
  TMatrixTSym<double> CovMatrix(NumberTags);
  for(std::size_t i = 0; i < NumberTags; i++) {
    CovMatrix(i, i) = m_Uncertainties[i]*m_Uncertainties[i];
  }
  // Set up multidimensional Gaussian containing predicted yields and normalized yields
  RooMultiVarGaussian Model("Model", "", m_NormalizedYields, m_PredictedYields, CovMatrix);
  // Set up dataset
  RooDataSet Data("Data", "", m_NormalizedYields);
  Data.add(m_NormalizedYields);
  // Perform fit
  auto Result = Model.fitTo(Data, RooFit::Save());
  Result->Print("V");
}

void FPlusFitter::AddMeasurement(const std::string &TagMode) {
  // Get raw double tag yields and efficiency
  double DT_RawYield = m_Settings[TagMode + "_DT_Yield"].getD(TagMode + "_DoubleTag_Yield");
  double DT_RawYield_err = m_Settings[TagMode + "_DT_Yield"].getD(TagMode + "_DoubleTag_Yield_err");
  double DT_Efficiency = m_Settings["DT_Efficiencies"].getD(TagMode + "_DoubleTagEfficiency");
  // Divide by efficiency
  double DT_Yield = DT_RawYield/DT_Efficiency;
  double DT_Yield_err = DT_RawYield_err/DT_Efficiency;
  // Get raw single tag yields and efficiency
  double ST_RawYield = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield");
  double ST_RawYield_err = m_Settings[TagMode + "_ST_Yield"].getD(TagMode + "_SingleTag_Yield_err");
  double ST_Efficiency = m_Settings["ST_Efficiencies"].getD(TagMode + "_SingleTagEfficiency");
  // Divide by efficiency
  double ST_Yield = ST_RawYield/ST_Efficiency;
  double ST_Yield_err = ST_RawYield_err/ST_Efficiency;
  // Create variable with a normalized yield and add to dataset
  auto NormalizedYield = Unique::create<RooRealVar*>((TagMode + "_Normalized_Yield").c_str(), "", DT_Yield/ST_Yield);
  m_NormalizedYields.add(*NormalizedYield);
  // Save statistical uncertainty
  m_Uncertainties.push_back(TMath::Sqrt(TMath::Power(DT_Yield_err/DT_Yield, 2) + TMath::Power(ST_Yield_err/ST_Yield, 2))*DT_Yield/ST_Yield);
}

void FPlusFitter::AddPrediction(const std::string &TagMode) {
  double KKpipi_BF = m_Settings["BranchingFractions"].getD("KKpipi");
  double FPlus_Tag = m_Settings["FPlus_TagModes"].getD(TagMode);
  auto PredictedYield = Unique::create<RooFormulaVar*>((TagMode + "_Normalized_Yield_Prediction").c_str(), Form("%f*(1 - (2*%f - 1)*(2*@0 - 1))", KKpipi_BF, FPlus_Tag), RooArgList(m_FPlus));
  m_PredictedYields.add(*PredictedYield);
}
