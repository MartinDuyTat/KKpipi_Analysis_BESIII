// Martin Duy Tat 31st March 2021
/**
 * PrepareTagTree is an application that takes in BESIII events selected in BOSS and applies initial cuts to a TTree, which is saved to a separate file
 * @param 1 Tag Mode that we want to study
 * @param 2 Type of tag mode, input "ST" for single tag, "DTSignal" for signal side of a double tag and "DTTag" for the tag side of a double tag
 * @param 3 Number of input ROOT files
 * @param 4 Filename of the ROOT files with BESIII data, without the final number and .root part
 * @param 5 Name of TTree
 * @param 6 Filename of the output ROOT file
 * @param 7 If this is an MC sample, specify the luminosity scale, otherwise it's assumed to be 1
 */

#include<iostream>
#include<cstdlib>
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"Utilities.h"
#include"InitialCuts.h"
#include"ApplyCuts.h"

int main(int argc, char *argv[]) {
  if(argc != 7 && argc != 8) {
    std::cout << "Need 6 or 7 input aguments\n";
    return 0;
  }
  std::string TagMode(argv[1]), TagType(argv[2]);
  if(TagType != "ST" && TagType != "DTSignal" && TagType != "DTTag") {
    std::cout << "Tag type " << argv[2] << " not recognized\n";
    return 0;
  }
  std::cout << "Sample preparation of " << TagMode << " tags of type " << TagType << "\n";
  std::cout << "Reading cuts...\n";
  TagType.erase(0, 2);
  InitialCuts Cuts(TagMode, TagType);
  ApplyCuts applyCuts(Cuts.GetInitialCuts());
  std::cout << "Cuts ready\n";
  std::cout << "Loading TChain...\n";
  TChain Chain;
  Utilities::LoadChain(&Chain, std::atoi(argv[3]), std::string(argv[4]), std::string(argv[5]));
  std::cout << "Applying cuts...\n";
  TFile OutputFile(argv[6], "RECREATE");
  TTree *OutputTree = applyCuts(&Chain);
  if(argc == 8) {
    OutputTree->SetWeight(1.0/std::atof(argv[7]));
  }
  OutputTree->SetDirectory(&OutputFile);
  OutputTree->Write();
  OutputFile.Close();
  std::cout << "Cuts applied and events saved to file\n";
  return 0;
}
