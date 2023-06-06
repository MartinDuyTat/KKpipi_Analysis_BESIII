// Martin Duy Tat 4th October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixT.h"
#include"TRandom.h"
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
  double YieldSmear = 0.0;
  if(settings.get("Systematics") == "STYield") {
    if(Tag != "KLpi0" && Tag != "KeNu") {
      const double STYield_err = ParsedSTYields.at(YieldName + "_err");
      YieldSmear = gRandom->Gaus(0.0, STYield_err);
    }
  } else if(settings.get("Systematics") == "KLpi0STYield") {
    if(Tag == "KLpi0") {
      const double STYield_err = ParsedSTYields.at(YieldName + "_err");
      YieldSmear = gRandom->Gaus(0.0, STYield_err);
    }
  } else if(settings.get("Systematics") == "KeNuSTYield") {
    if(Tag == "KeNu") {
      const double STYield_err = ParsedSTYields.at(YieldName + "_err");
      YieldSmear = gRandom->Gaus(0.0, STYield_err);
    }
  } else if(settings.get("Systematics") == "PeakingBackgrounds") {
    const auto SystematicsYieldName = YieldName + "_PeakingBackgrounds_syst_err";
    if(ParsedSTYields.find(SystematicsYieldName) != ParsedSTYields.end()) {
      const double STPeakingBackground_err =
	ParsedSTYields.at(SystematicsYieldName);
      YieldSmear = gRandom->Gaus(0.0, STPeakingBackground_err);
      std::cout << YieldSmear << "\n";
    }
  }
  const std::string EffName = Tag + "_SingleTagEfficiency";
  double Efficiency = settings["ST_Efficiency"].getD(EffName);
  if(settings.get("Systematics") == "Efficiency") {
    const double Eff_err = settings["ST_Efficiency"].getD(EffName + "_err");
    Efficiency += gRandom->Gaus(0.0, Eff_err);
  }
  return (ParsedSTYields.at(YieldName) + YieldSmear)/Efficiency;
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
  if(settings.get("Systematics") == "Efficiency") {
    TMatrixT<double> *EffMatrix_err = nullptr;
    EffMatrixFile.GetObject((EffMatrixName + "_err").c_str(), EffMatrix_err);
    for(int i = 0; i < EffMatrix->GetNrows(); i++) {
      for(int j = 0; j < EffMatrix->GetNcols(); j++) {
	if((*EffMatrix_err)(i, j) != 0.0) {
	  (*EffMatrix)(i, j) += gRandom->Gaus(0.0, (*EffMatrix_err)(i, j));
	  if((*EffMatrix)(i, j) < 0.0) {
	    (*EffMatrix)(i, j) = 0.0;
	  }
	}
      }
    }
  }
  EffMatrixFile.Close();
  return *EffMatrix;
}
