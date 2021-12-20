// Martin Duy Tat 4th April 2021
/**
 * FitSingleTagMBC is an application that determines the single tag yield using a fit of the beam constrained mass \f$m_\text{BC}\f$ distribution
 * The signal shape is taken from an exclusive signal MC sample, convolved with a Gaussian, and the background is modelled with an Argus PDF
 */

#include<iostream>
#include<string>
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"SingleTagYield.h"
#include"Utilities.h"
#include"Settings.h"

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "Single tag yield fit\n";
  std::cout << "Loading ROOT files...\n";
  std::string TreeName = settings.get("TreeName");
  TChain Chain(TreeName.c_str());
  std::string Mode = settings.get("Mode");
  std::string DataFilename = Utilities::ReplaceString(settings["Datasets_WithDeltaECuts"].get("Dskim"), "TAG", Mode);
  Chain.Add(DataFilename.c_str());
  std::string SignalMCFilename = settings["Datasets_WithDeltaECuts"].get("SignalMC_ST");
  SignalMCFilename = Utilities::ReplaceString(SignalMCFilename, "TAG", Mode);
  TFile MCFile(SignalMCFilename.c_str(), "READ");
  TTree *MCTree = nullptr;
  MCFile.GetObject(TreeName.c_str(), MCTree);
  MCTree->SetBranchStatus("*", 0);
  MCTree->SetBranchStatus("MBC", 1);
  std::cout << "Trees ready\n";
  std::cout << "Setting up fit model...\n";
  TTree *ClonedMCTree = nullptr;
  if(settings.getI("Events_in_MC") < 0 || settings.getI("Events_in_MC") > MCTree->GetEntries()) {
    ClonedMCTree = MCTree;
  } else {
    ClonedMCTree = MCTree->CloneTree(settings.getI("Events_in_MC"));
  }
  SingleTagYield singleTagYield(&Chain, ClonedMCTree, settings);
  std::cout << "Ready to fit\n";
  std::cout << "Fitting single tag yield of " << Mode << "...\n";
  singleTagYield.FitYield();
  std::cout << "Fit completed" << "\n";
  return 0;
}
