// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixT.h"
#include"BinnedDTYieldPrediction.h"
#include"Settings.h"
#include"Utilities.h"

BinnedDTYieldPrediction::BinnedDTYieldPrediction(const std::string &Tag,
						 const Settings &settings):
  m_SingleTagYield(GetSTYield(Tag, settings)),
  m_EfficiencyMatrix(GetEffMatrix(Tag, settings, "EffMatrix")),
  m_EfficiencyMatrix_CPEven(GetEffMatrix(Tag, settings, "EffMatrix_CPEven")),
  m_EfficiencyMatrix_CPOdd(GetEffMatrix(Tag, settings, "EffMatrix_CPOdd")),
  m_EfficiencyMatrix_K0pipi(GetEffMatrix(Tag, settings, "EffMatrix_K0pipi")),
  m_y(settings["Mixing"].getD("y")) {
}

double BinnedDTYieldPrediction::GetSTYield(const std::string &Tag,
					   const Settings &settings) const {
  const std::string STYieldFilename =
    Utilities::ReplaceString(settings.get("ST_Yield"), "TAG", Tag);
  const auto ParsedSTYields = Utilities::ParseFile(STYieldFilename);
  const std::string YieldName(Tag + "_SingleTag_Yield");
  const std::string EffName = Tag + "_SingleTagEfficiency";
  const double Efficiency = settings["ST_Efficiency"].getD(EffName);
  return ParsedSTYields.at(YieldName)/Efficiency;
}

TMatrixT<double> BinnedDTYieldPrediction::GetEffMatrix(
  const std::string &Tag,
  const Settings &settings,
  const std::string &EffMatrixName) const {
  const std::string Filename =
    Utilities::ReplaceString(settings.get("EffMatrix"), "TAG", Tag);
  TFile EffMatrixFile(Filename.c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  EffMatrixFile.GetObject(EffMatrixName.c_str(), EffMatrix);
  /*if(Smearing && m_Settings.get("Systematics") == "Efficiency") {
    TMatrixT<double> *EffMatrix_err = nullptr;
    EffMatrixFile.GetObject("EffMatrix_err", EffMatrix_err);
    for(int i = 0; i < EffMatrix->GetNrows(); i++) {
      for(int j = 0; j < EffMatrix->GetNcols(); j++) {
	(*EffMatrix)(i, j) += gRandom->Gaus(0.0, (*EffMatrix_err)(i, j));
      }
    }
  }*/
  EffMatrixFile.Close();
  //EffMatrix->Invert();
  return *EffMatrix;
}
