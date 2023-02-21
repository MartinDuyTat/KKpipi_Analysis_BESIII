// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixT.h"
#include"BinnedDTYieldPrediction.h"
#include"Settings.h"
#include"FlavourTags/KiCombiner.h"

BinnedDTYieldPrediction::BinnedDTYieldPrediction(const std::string &Tag,
						 const KiCombiner *Ki,
						 const Settings &settings):
  m_SingleTagYield(GetSTYield(Tag, settings)),
  m_EfficiencyMatrix_CPEven(GetEffMatrix(Tag, settings, +1)),
  m_EfficiencyMatrix_CPOdd(GetEffMatrix(Tag, settings, -1)),
  m_Ki(Ki) {
}

double BinnedDTYieldPrediction::GetSTYield(const std::string &Tag,
					   const Settings &settings) const {
  const std::string SettingsName(Tag + "_ST_Yield");
  const std::string YieldName(Tag + "_SingleTag_Yield");
  const std::string EffName = Tag + "_SingleTagEfficiency";
  const double Efficiency = settings["ST_Efficiency"].getD(EffName);
  return settings[SettingsName].getD(YieldName)/Efficiency;
}

TMatrixT<double> BinnedDTYieldPrediction::GetEffMatrix(
  const std::string &Tag, const Settings &settings, int CPEvenOdd) const {
  TFile EffMatrixFile(settings.get(Tag + "_EfficiencyMatrix").c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  const std::string EffMatrixName = CPEvenOdd > 0 ?
                                    "EffMatrix_CPEven" :
                                    "EffMatrix_CPOdd";
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
