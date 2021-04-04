// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include<iostream>
#include"TCut.h"
#include"InitialCuts.h"

InitialCuts::InitialCuts(const std::string &TagMode, const std::string &TagType): m_TagMode(TagMode), m_TagType(TagType) {
}

TCut InitialCuts::GetInitialCuts() const {
  if(m_TagType == "ST") {
    return GetTagCuts(m_TagMode);
  } else if(m_TagType == "DT") {
    return GetTagCuts("KKpipi", "Signal") && GetTagCuts(m_TagMode, "Tag");
  } else {
    return TCut();
  }
}

TCut InitialCuts::GetTagCuts(const std::string &TagMode, const std::string &TagSide) const {
  std::ifstream CutFile(std::string(INITIAL_CUTS_DIR) + TagMode + ".cut");
  TCut Cuts;
  if(CutFile.is_open()) {
    std::string line;
    while(std::getline(CutFile, line)) {
      if(line != "No cuts") {
	line = ReplaceTagMode(line, TagSide);
	Cuts = Cuts && TCut(line.c_str());
      } else {
	return TCut();
      }
    }
  } else {
    std::cout << "Could not find file " << std::string(INITIAL_CUTS_DIR) + TagMode + ".cut\n";
  }
  return Cuts;
}

std::string InitialCuts::ReplaceTagMode(std::string Cuts, const std::string &TagSide) const {
  std::string ReplaceThis("SignalTag");
  while(true) {
    std::string::size_type Position = Cuts.find(ReplaceThis);
    if(Position == std::string::npos) {
      break;
    } else {
      Cuts.replace(Position, ReplaceThis.length(), TagSide);
    }
  }
  return Cuts;
}
