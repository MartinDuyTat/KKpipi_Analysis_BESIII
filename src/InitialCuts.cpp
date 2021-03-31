// Martin Duy Tat 31st March 2021

#include<string>
#include<fstream>
#include<iostream>
#include"TCut.h"
#include"InitialCuts.h"

InitialCuts::InitialCuts(const std::string &TagMode, const std::string &TagType): m_TagMode(TagMode), m_TagType(TagType) {
}

TCut InitialCuts::GetInitialCuts() const {
  std::ifstream CutFile(INITIAL_CUTS_DIR + "/InitialCuts/" + m_TagMode + ".cut");
  TCut Cuts;
  if(CutFile.is_open()) {
    std::string line;
    while(std::getline(CutFile, line)) {
      if(line != "No cuts") {
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
