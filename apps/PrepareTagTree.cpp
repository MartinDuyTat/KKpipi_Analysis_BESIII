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
#include"Settings.h"

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::string Mode = settings.get("Mode");
  std::string TagType = settings.get("TagType");
  bool IncludeDeltaECuts = settings.getB("Include_DeltaE_Cuts");
  bool TruthMatch = settings.getB("TruthMatch");
  std::string TreeName = settings.get("TreeName");
  int DataSetType = settings.getI("DataSetType");
  double LuminosityScale = settings.getD("LuminosityScale");
  std::string InputFiles = settings.get("InputFiles");
  std::string OutputFilename = settings.get("OutputFilename");
  std::cout << "Sample prepration of " << TagType << " " << Mode << " mode\n";
  std::cout << "DataSetType: " << DataSetType << "\n";
  std::string DataMC = DataSetType == 0 ? "Data" : "MC";
  TCut Cuts = Utilities::LoadCuts(Mode, IncludeDeltaECuts, TruthMatch, TagType, DataMC);
  // This cut removes empty NTuples
  Cuts = Cuts && TCut("!(Run == 0 && Event == 0)");
  std::cout << "Cuts ready, will apply the following cuts:\n" << Cuts.GetTitle() << "\n";
  ApplyCuts applyCuts(Cuts);
  std::cout << "Cuts ready\n";
  std::cout << "Loading TChain...\n";
  TChain Chain;
  Utilities::LoadChain(&Chain, InputFiles, TreeName);
  std::cout << "Applying cuts...\n";
  TFile OutputFile(OutputFilename.c_str(), "RECREATE");
  TTree *OutputTree = applyCuts(&Chain, DataSetType, LuminosityScale);
  OutputTree->SetDirectory(0);
  OutputTree->Write();
  OutputFile.Close();
  //delete OutputTree;
  std::cout << "Cuts applied and events saved to file " << OutputFilename << "\n";
  return 0;
}
