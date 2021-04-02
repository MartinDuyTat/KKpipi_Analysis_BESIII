// Martin Duy Tat 31st March 2021

#include<iostream>
#include<fstream>
#include<string>
#include"TChain.h"
#include"Utilities.h"

namespace Utilities {
  void LoadChain(TChain *Chain, int NumberFiles, const std::string &Filename, const std::string &TreeName) {
    std::cout << "Initializing TChain with files...\n";
    if(TreeName != "") {
      Chain->SetName(TreeName.c_str());
    }
    if(NumberFiles == 0) {
      std::cout << "Need more than 0 input files...\n";
      return;
    }
    for(int i = 0; i < NumberFiles; i++) {
      Chain->Add((Filename + std::to_string(i) + ".root").c_str());
    }
    std::cout << "ROOT files added to TChain\n";
    return;
  }

  void LoadChain(TChain *Chain, const std::string &Filename, const std::string &TreeName) {
    std::cout << "Initializing TChain with files...\n";
    if(TreeName != "") {
      Chain->SetName(TreeName.c_str());
    }
    std::ifstream Infile(Filename);
    std::string line;
    while(std::getline(Infile, line)) {
      Chain->Add(line.c_str());
    }
    Infile.close();
    std::cout << "ROOT files added to TChain\n";
    return;
  }
}
