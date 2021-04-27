// Martin Duy Tat 27th April 2021

#include<fstream>
#include<iostream>
#include"CutsFromFile.h"

std::string CutsFromFile::ReplaceTagMode(std::string Cuts, const std::string &TagSide) const {
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

TCut CutsFromFile::GetTagCuts(const std::string &Filename, const std::string &TagSide) const {
  std::ifstream CutFile(Filename);
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
    std::cout << "Could not find file " << Filename << "\n";
  }
  return Cuts;
}

CutsFromFile::CutsFromFile(const std::string &Filename, const std::string &TagSide ) {
  m_Cuts = GetTagCuts(Filename, TagSide);
}

TCut CutsFromFile::GetCuts() const {
  return m_Cuts;
}
