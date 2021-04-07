// Martin Duy Tat 7th April 2021

#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cctype>
#include<tuple>
#include<utility>
#include"TFile.h"
#include"TTree.h"
#include"TopoAnaReader.h"

TopoAnaReader::TopoAnaReader(const std::string &Filename) {
  std::ifstream Infile(Filename);
  std::string line;
  m_OtherComponent.first = "Other";
  while(std::getline(Infile, line)) {
    std::string::size_type DecayStartPosition = line.find('(');
    if(DecayStartPosition != std::string::npos) {
      std::string ComponentLabel = line.substr(0, DecayStartPosition);
      std::string DecayDescriptor = line.substr(DecayStartPosition + 1);
      std::string::size_type DecayEndPosition = DecayDescriptor.find(')');
      if(DecayEndPosition == std::string::npos) {
	continue;
      }
      DecayDescriptor = DecayDescriptor.substr(0, DecayEndPosition);
      m_DecayComponents.push_back(std::make_tuple(ComponentLabel, DecayDescriptor, std::vector<int>()));
    } else {
      std::string::iterator EndPosition = std::remove(line.begin(), line.end(), ' ');
      line.erase(EndPosition, line.end());
      m_OtherComponent.first = line;
    }
  }
  Infile.close();
}

void TopoAnaReader::AnalyzeComponents(const std::string &Filename) {
  std::ifstream Infile(Filename);
  std::string line;
  while(std::getline(Infile, line)) {
    std::stringstream ss(line);
    std::string word;
    ss >> word;
    if(word != "rowNo:") {
      continue;
    }
    ss >> word >> word >> word;
    int DecayTopology = std::stoi(word);
    std::string D0Decay, D0barDecay;
    std::getline(Infile, line);
    std::getline(Infile, D0Decay);
    std::getline(Infile, D0barDecay);
    bool FoundDecay = false;
    for(auto DecayComponent_iter = m_DecayComponents.begin(); DecayComponent_iter != m_DecayComponents.end(); DecayComponent_iter++) {
      std::string DecayDescriptor = std::get<1>(*DecayComponent_iter);
      std::string::size_type FoundD0Decay = D0Decay.find(DecayDescriptor + "  &");
      std::string::size_type FoundD0barDecay = D0barDecay.find(DecayDescriptor + "  &");
      if(FoundD0Decay != std::string::npos || FoundD0barDecay != std::string::npos) {
	std::get<2>(*DecayComponent_iter).push_back(DecayTopology);
	FoundDecay = true;
	break;
      }
    }
    if(!FoundDecay) {
      m_OtherComponent.second.push_back(DecayTopology);
    }
  }
}

void TopoAnaReader::SaveAllComponentCuts() const {
  for(auto Component_iter = m_DecayComponents.begin(); Component_iter != m_DecayComponents.end(); Component_iter++) {
    SaveComponentCuts(std::get<0>(*Component_iter), std::get<2>(*Component_iter));
  }
  SaveComponentCuts(m_OtherComponent.first, m_OtherComponent.second);
}

void TopoAnaReader::SaveComponentCuts(const std::string &Label, const std::vector<int> &DecayTopologies) const {
  std::ofstream Outfile(Label + ".cut");
  for(auto iter = DecayTopologies.begin(); iter != DecayTopologies.end(); iter++) {
    if(iter != DecayTopologies.begin()) {
      Outfile << " || ";
    }
    Outfile << "iDcyTr == " << *iter;
  }
}

void TopoAnaReader::SaveTree(TTree *InTree, const std::string &Label, const std::vector<int> &DecayTopologies) const {
  TFile OutputFile((Label + std::string(".root")).c_str(), "RECREATE");
  TTree *Tree = InTree->CloneTree(0);
  int iDcyTr;
  InTree->SetBranchAddress("iDcyTr", &iDcyTr);
  for(int i = 0; i < InTree->GetEntries(); i++) {
    InTree->GetEntry(i);
    if(std::find(DecayTopologies.begin(), DecayTopologies.end(), iDcyTr) != DecayTopologies.end()) {
      Tree->Fill();
    }
  }
  Tree->Write();
  OutputFile.Close();
}
    

void TopoAnaReader::SaveAllTrees(TTree *InTree) const {
  for(auto iter = m_DecayComponents.begin(); iter != m_DecayComponents.end(); iter++) {
    if(std::get<2>(*iter).size() != 0) {
      SaveTree(InTree, std::get<0>(*iter), std::get<2>(*iter));
    }
  }
  if(m_OtherComponent.second.size() != 0) {
    SaveTree(InTree, m_OtherComponent.first, m_OtherComponent.second);
  }
}
