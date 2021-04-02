// Martin Duy Tat 2nd April 2021

#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<algorithm>
#include"TCut.h"
#include"DeltaECut.h"

DeltaECut::DeltaECut(const std::string &TagMode, const std::string &TagType): m_TagMode(TagMode), m_TagType(TagType) {
}

TCut DeltaECut::GetDeltaECut() const {
  std::ifstream CutFile(std::string(DELTAE_CUTS_DIR) + m_TagMode + ".cut");
  if(CutFile.is_open()) {
    std::string line;
    std::getline(CutFile, line);
    double Lower, Upper;
    std::stringstream ss(line);
    ss >> Lower >> Upper;
    if(Lower > Upper) {
      std::swap(Lower, Upper);
    }
    std::string Cut = std::string(m_TagType + "DeltaE > " + std::to_string(Lower) + " && " + m_TagType + "DeltaE < " + std::to_string(Upper));
    return TCut(Cut.c_str());
  } else {
    std::cout << "Could not find file " << std::string(INITIAL_CUTS_DIR) + m_TagMode + ".cut\n";
    return TCut();
  }
}
