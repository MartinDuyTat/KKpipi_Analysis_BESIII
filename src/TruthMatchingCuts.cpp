// Martin Duy Tat 4th April 2021

#include<string>
#include"TCut.h"
#include"TruthMatchingCuts.h"
#include"CutsFromFile.h"

TruthMatchingCuts::TruthMatchingCuts(const std::string &TagMode, const std::string &TagType): m_TagMode(TagMode), m_TagType(TagType) {
}

TCut TruthMatchingCuts::GetModeSpecificCuts() const {
  if(m_TagType == "ST") {
    CutsFromFile cutsFromFile(std::string(TRUTH_MATCHING_CUTS_DIR) + m_TagMode + ".cut", "");
    return cutsFromFile.GetCuts();
  } else if (m_TagType == "DT") {
    CutsFromFile cutsFromFileSignal(std::string(TRUTH_MATCHING_CUTS_DIR) + "KKpipi" + ".cut", "Signal");
    CutsFromFile cutsFromFileTag(std::string(TRUTH_MATCHING_CUTS_DIR) + m_TagMode + ".cut", "Tag");
    return cutsFromFileSignal.GetCuts() && cutsFromFileTag.GetCuts();
  } else {
    return TCut();
  }
}

TCut TruthMatchingCuts::GetTruthMatchingCuts() const {
  if(m_TagType == "ST") {
    return TCut("IsSameDMother == 1 && PIDTrue == 1") && GetModeSpecificCuts();
  } else if(m_TagType == "DT") {
    return TCut("SignalIsSameDMother == 1 && SignalPIDTrue == 1 && TagIsSameDMother == 1 && TagPIDTrue == 1") && GetModeSpecificCuts();
  } else {
    return TCut();
  }
}
