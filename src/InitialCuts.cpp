// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include<iostream>
#include"TCut.h"
#include"InitialCuts.h"

InitialCuts::InitialCuts(const std::string &TagMode, const std::string &TagType): m_TagMode(TagMode), m_TagType(TagType) {
}

TCut InitialCuts::GetInitialCuts() const {
  std::ifstream CutFile(INITIAL_CUTS_DIR + m_TagMode + ".cut");
  TCut Cuts;
  if(CutFile.is_open()) {
    std::string line;
    while(std::getline(CutFile, line)) {
      if(line != "No cuts") {
	line = ReplaceTagMode(line);
	Cuts = Cuts && TCut(line.c_str());
      } else {
	return TCut();
      }
    }
  } else {
    std::cout << "Could not find file " << INITIAL_CUTS_DIR + "/InitialCuts/" + m_TagMode + ".cut\n";
  }
  return Cuts;
}

std::string InitialCuts::ReplaceTagMode(std::string Cuts) {
  std::string ReplaceThis("SignalTag");
  while(true) {
    std::string::size_type Position = Cuts.find(ReplaceThis);
    if(Position == std::string::npos) {
      break;
    } else {
      Cuts.replace(Position, ReplaceThis.length(), m_TagType);
    }
  }
  return Cuts;
}
