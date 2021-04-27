// Martin Duy Tat 7th April 2021
/**
 * AnalyzeTopoAna takes in a list of particular decay components in a TopoAna file and finds all the decay topology numbers that match those decay components, the rest are classified as "Other" decay components
 * @param 1 Filename of TopoAna txt output file
 * @param 2 Filename with decay components
 * @param 3 Filename of original ROOT file
 * @param 4 Name of TTree
 */

#include<iostream>
#include"TFile.h"
#include"TChain.h"
#include"TopoAnaReader.h"

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "Need 4 input arguments\n";
    return 0;
  }
  std::cout << "Loading the decay components...\n";
  TopoAnaReader Reader{std::string(argv[2])};
  std::cout << "Decay components ready\n";
  std::cout << "Going through TopoAna output...\n";
  Reader.AnalyzeComponents(std::string(argv[1]));
  std::cout << "TopoAna analysis complete\n";
  std::cout << "Saving cut files...\n";
  Reader.SaveAllComponentCuts();
  std::cout << "Cuts safely stored in files\n";
  std::cout << "Saving TTrees...\n";
  TChain Chain(argv[4]);
  Chain.Add(argv[3]);
  Reader.SaveAllTrees(&Chain);
  std::cout << "All trees saved\n";
  return 0;
}
