// Martin Duy Tat 31st March 2021
/**
 * PrepareTagTree is an application that takes in BESIII events selected in BOSS and applies initial cuts to a TTree, which is saved to a separate file
 * @param 1 Tag Mode that we want to study
 * @param 2 Type of tag mode, input "ST" for single tag, "DT" for double tag
 * @param 3 Type "DeltaECuts" to use standard initial cuts plus \f$\Delta E\f$ cuts, type "NoDeltaECuts" to use standard inital cuts only, and "TruthMatchingCuts" to truth match the sample
 * @param 4 Number of input ROOT files
 * @param 5 Filename of the ROOT files with BESIII data, without the final number and .root part
 * @param 6 Name of TTree
 * @param 7 Filename of the output ROOT file
 * @param 8 If this is an MC sample, specify the luminosity scale, otherwise it's assumed to be 1
 */

#include<iostream>
#include<cstdlib>
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"Utilities.h"
#include"ApplyCuts.h"

int main(int argc, char *argv[]) {
  if(argc != 8 && argc != 9) {
    std::cout << "Need 7 or 8 input aguments\n";
    return 0;
  }
  std::string TagMode(argv[1]), TagType(argv[2]), CutType(argv[3]);
  if(TagType != "ST" && TagType != "DT") {
    std::cout << "Tag type " << argv[2] << " not recognized\n";
    return 0;
  }
  std::cout << "Sample preparation of " << TagMode << " tags of type " << TagType << "\n";
  std::cout << "Reading cuts...\n";
  TCut Cuts = Utilities::LoadCuts(CutType, TagMode, TagType);
  std::cout << "Cuts ready, will apply the following cuts:\n" << Cuts.GetTitle() << "\n";
  ApplyCuts applyCuts(Cuts);
  std::cout << "Cuts ready\n";
  std::cout << "Loading TChain...\n";
  TChain Chain;
  Utilities::LoadChain(&Chain, std::atoi(argv[4]), std::string(argv[5]), std::string(argv[6]));
  std::cout << "Applying cuts...\n";
  TFile OutputFile(argv[7], "RECREATE");
  TTree *OutputTree = applyCuts(&Chain);
  if(argc == 9) {
    double Weight = 1.0/std::atof(argv[8]);
    std::cout << "Weighting the events by " << Weight << " to account for luminosity scale\n";
    OutputTree->SetWeight(Weight);
  }
  OutputTree->SetDirectory(&OutputFile);
  OutputTree->Write();
  OutputFile.Close();
  std::cout << "Cuts applied and events saved to file\n";
  return 0;
}
