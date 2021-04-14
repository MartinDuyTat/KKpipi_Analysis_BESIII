// Martin Duy Tat 31st March 2021
/**
 * PrepareTagTree is an application that takes in BESIII events selected in BOSS and applies initial cuts to a TTree, which is saved to a separate file
 * Multiple datasets can be processed, using a parsed options file, the application will look for a line with "Name" and parse all the options from the subsequent lines
 * Each option must be on a separate line, in the following format and order:
 * Line 1: Name <Name of dataset> (not used for anything but printing messages)
 * Line 2: Number of ROOT files
 * Line 3: Filename of the ROOT files with BESIII data, without the final number and .root part
 * Line 4: Tag Mode that we want to study
 * Line 5: Type of tag mode, input "ST" for single tag, "DT" for double tag
 * Line 6: Type "DeltaECuts" to use standard initial cuts plus \f$\Delta E\f$ cuts, type "NoDeltaECuts" to use standard inital cuts only, and "TruthMatchingCuts" to truth match the sample
 * Line 7: Name of TTree
 * Line 8: Write a number between \f$0\f$ and \f$9\f$, where each number represents, in this order: Data, \f$D^0\bar{D^0}\f$ MC, \f$D^+D^-\f$ MC, $q\bar{q}$ MC, \f$\psi'\gamma\f$ MC, \f$J/\psi\gamma\f$ MC, \f$\tau\tau\f$ MC, \f$\mu\mu\f$ MC, \f$ee\f$ MC, non-\f$D\bar{D}\f$ MC
 * Line 9: Luminosity scale of MC, for data input 1
 * Line 10: Filename of the output ROOT file
 * @param 1 Filename of options file
 */

#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"Utilities.h"
#include"ApplyCuts.h"

struct TagTreeSettings {
  std::string Name;
  int NumberFiles;
  std::string InputFilename;
  std::string TagMode;
  std::string TagType;
  std::string CutType;
  std::string TreeName;
  int DataSetType;
  double LuminosityScale;
  std::string OutputFilename;
};

void PrepareTagTree(TagTreeSettings Settings);

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input agument\n";
    return 0;
  }
  std::ifstream Infile(argv[1]);
  if(!Infile.is_open()) {
    std::cout << "Unable to open options file\n";
    return 0;
  }
  while(Infile.peek() != EOF) {
    std::string First;
    Infile >> First;
    if(First == "Name") {
      TagTreeSettings Settings;
      Infile >> Settings.Name;
      Infile >> Settings.NumberFiles;
      Infile >> Settings.InputFilename;
      Infile >> Settings.TagMode;
      Infile >> Settings.TagType;
      Infile >> Settings.CutType;
      Infile >> Settings.TreeName;
      Infile >> Settings.DataSetType;
      Infile >> Settings.LuminosityScale;
      Infile >> Settings.OutputFilename;
      PrepareTagTree(Settings);
    }
  }
  return 0;
}

void PrepareTagTree(TagTreeSettings Settings) {
  std::cout << "Sample prepartion of " << Settings.Name << " datset\n";
  if(Settings.TagType != "ST" && Settings.TagType != "DT") {
    std::cout << "Tag type " << Settings.TagType << " not recognized\n";
    return;
  }
  std::cout << "Tag mode: " << Settings.TagMode << "\nTag type: " << Settings.TagType << "\n";
  std::cout << "Reading cuts...\n";
  std::string DataMC;
  if(Settings.DataSetType == 0) {
    DataMC = "Data";
  } else if(Settings.DataSetType > 0 && Settings.DataSetType < 10) {
    DataMC = "MC";
  } else {
    std::cout << "DataSetType " << Settings.DataSetType << " not recognized\n";
    return;
  }
  TCut Cuts = Utilities::LoadCuts(Settings.CutType, Settings.TagMode, Settings.TagType, DataMC);
  // This cut removes empty NTuples
  //Cuts = Cuts && TCut("!(Run == 0 && Event == 0)");
  Cuts = Cuts && TCut("MBC != 0");
  std::cout << "Cuts ready, will apply the following cuts:\n" << Cuts.GetTitle() << "\n";
  ApplyCuts applyCuts(Cuts);
  std::cout << "Cuts ready\n";
  std::cout << "Loading TChain...\n";
  TChain Chain;
  Utilities::LoadChain(&Chain, Settings.NumberFiles, Settings.InputFilename, Settings.TreeName);
  std::cout << "Applying cuts...\n";
  TFile OutputFile(Settings.OutputFilename.c_str(), "RECREATE");
  TTree *OutputTree = applyCuts(&Chain, Settings.DataSetType, Settings.LuminosityScale);
  OutputTree->SetDirectory(&OutputFile);
  OutputTree->Write();
  OutputFile.Close();
  std::cout << "Cuts applied and events saved to file " << Settings.OutputFilename << "\n";
}
