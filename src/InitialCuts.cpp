// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include<iostream>
#include"TCut.h"
#include"InitialCuts.h"
#include"CutsFromFile.h"

InitialCuts::InitialCuts(const std::string &TagMode, const std::string &TagType, bool KKpipiPartReco): m_TagMode(TagMode), m_TagType(TagType), m_KKpipiPartReco(KKpipiPartReco) {
}

TCut InitialCuts::GetInitialCuts() const {
  if(m_TagType == "ST") {
    return GetTagCuts(m_TagMode);
  } else if(m_TagType == "DT") {
    if(!m_KKpipiPartReco) {
      return GetTagCuts("KKpipi", "Signal") && GetTagCuts(m_TagMode, "Tag");
    } else {
      return GetTagCuts("KKpipiPartReco", "Signal") && GetTagCuts(m_TagMode, "Tag");
    }
  } else {
    return TCut();
  }
}

TCut InitialCuts::GetTagCuts(const std::string &TagMode, const std::string &TagSide) const {
  CutsFromFile cutsFromFile(std::string(INITIAL_CUTS_DIR) + TagMode + ".cut", TagSide);
  return cutsFromFile.GetCuts();
}
