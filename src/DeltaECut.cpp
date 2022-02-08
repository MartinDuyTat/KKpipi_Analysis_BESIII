// Martin Duy Tat 2nd April 2021

#include<string>
#include<iostream>
#include<sstream>
#include<algorithm>
#include<stdexcept>
#include"TCut.h"
#include"DeltaECut.h"
#include"Settings.h"

DeltaECut::DeltaECut(const std::string &TagMode,
		     const std::string &TagType,
		     const std::string &DataMC,
		     bool KKpipiPartReco): m_TagMode(TagMode),
					   m_TagType(TagType),
					   m_DataMC(DataMC),
					   m_KKpipiPartReco(KKpipiPartReco) {
  if(m_TagType == "ST") {
    if(m_TagMode.find("KL") != std::string::npos || m_TagMode == "KeNu" || m_KKpipiPartReco) {
      throw std::invalid_argument("Cannot have partially reconstructed single tag");
    }
  }
}

TCut DeltaECut::GetDeltaECut() const {
  if(m_TagType == "ST") {
    return GetDeltaECutFromFile(m_TagMode);
  } else if(m_TagType == "DT") {
    if(m_TagMode.find("KL") != std::string::npos || m_TagMode == "KeNu") {
      return GetDeltaECutFromFile("KKpipi", "Signal");
    } else if(m_KKpipiPartReco) {
      return GetDeltaECutFromFile(m_TagMode, "Tag");
    } else{
      return GetDeltaECutFromFile("KKpipi", "Signal") && GetDeltaECutFromFile(m_TagMode, "Tag");
    }
  } else {
    return TCut();
  }
}

TCut DeltaECut::GetDeltaECutFromFile(const std::string &TagMode, const std::string &TagSide) const {
  std::string CutFilename(std::string(DELTAE_CUTS_DIR) + "DeltaECuts_" + m_DataMC + ".cut");
  Settings DeltaESettings("DeltaESettings", CutFilename);
  double Lower = DeltaESettings.getD(TagMode + "_DeltaE_LowerCut");
  double Upper = DeltaESettings.getD(TagMode + "_DeltaE_UpperCut");
  std::string Cut = std::string(TagSide + "DeltaE > " + std::to_string(Lower) + " && " + TagSide + "DeltaE < " + std::to_string(Upper));
  return TCut(Cut.c_str());
}
