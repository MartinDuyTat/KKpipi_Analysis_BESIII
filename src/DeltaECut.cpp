// Martin Duy Tat 2nd April 2021

#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<algorithm>
#include"TCut.h"
#include"DeltaECut.h"

DeltaECut::DeltaECut(const std::string &TagMode, const std::string &TagType, const std::string &DataMC): m_TagMode(TagMode), m_TagType(TagType), m_DataMC(DataMC) {
}

TCut DeltaECut::GetDeltaECut() const {
  if(m_TagType == "ST") {
    return GetDeltaECutFromFile(m_TagMode);
  } else if(m_TagType == "DT") {
    return GetDeltaECutFromFile("KKpipi", "Signal") && GetDeltaECutFromFile(m_TagMode, "Tag");
  } else {
    return TCut();
  }
}

TCut DeltaECut::GetDeltaECutFromFile(const std::string &TagMode, const std::string &TagSide) const {
  std::ifstream CutFile(std::string(DELTAE_CUTS_DIR) + TagMode + m_DataMC + ".cut");
  if(CutFile.is_open()) {
    std::string line;
    std::getline(CutFile, line);
    double Lower, Upper;
    std::stringstream ss(line);
    ss >> Lower >> Upper;
    if(Lower > Upper) {
      std::swap(Lower, Upper);
    }
    std::string Cut = std::string(TagSide + "DeltaE > " + std::to_string(Lower) + " && " + TagSide + "DeltaE < " + std::to_string(Upper));
    return TCut(Cut.c_str());
  } else {
    std::cout << "Could not find file " << std::string(DELTAE_CUTS_DIR) + TagMode + m_DataMC + ".cut\n";
    return TCut();
  }
}
