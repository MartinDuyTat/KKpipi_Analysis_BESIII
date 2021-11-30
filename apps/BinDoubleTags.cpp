// Martin Duy Tat 26th November 2021
/**
 * BinDoubleTags is an application that determines which phase space bin a double tag belong to
 * Both the signal and tag side are analyzed and a combined phase space bin is saved to the ROOT file
 */

#include<iostream>
#include<memory>
#include<string>
#include<utility>
#include"TChain.h"
#include"TTree.h"
#include"TFile.h"
#include"Utilities.h"
#include"Settings.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"

int main(int argc, char *argv[]) {
  std::cout << "Binning double tag events\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::string TreeName = settings.get("TreeName");
  std::cout << "Loading double tag events and setting up output TTree...\n";
  TChain InputChain(TreeName.c_str());
  InputChain.Add(settings.get("InputFilename").c_str());
  std::unique_ptr<KKpipi_PhaseSpace> PhaseSpace = Utilities::GetPhaseSpaceBinning(settings, &InputChain);
  std::string OutputFilename = settings.get("OutputFilename");
  std::string SignalBin_Name = settings.get("SignalBin_variable");
  std::string TagBin_Name = settings.get("TagBin_variable");
  int SignalBin, TagBin, SignalBin_true, TagBin_true;
  TFile OutputFile(OutputFilename.c_str(), "RECREATE");
  TTree *OutputTree = InputChain.CloneTree(0);
  if(settings.getB("Bin_reconstructed")) {
    OutputTree->Branch(SignalBin_Name.c_str(), &SignalBin);
    OutputTree->Branch(TagBin_Name.c_str(), &TagBin);
  }
  if(settings.getB("Bin_truth")) {
    OutputTree->Branch((SignalBin_Name + "_true").c_str(), &SignalBin_true);
    OutputTree->Branch((TagBin_Name + "_true").c_str(), &TagBin_true);
  }
  int EventsOutsidePhaseSpace = 0;
  int EventsOutsidePhaseSpace_true = 0;
  std::cout << "Ready to bin phase space\n";
  for(int i = 0; i < InputChain.GetEntries(); i++) {
    InputChain.GetEntry(i);
    if(settings.getB("Bin_reconstructed")) {
      std::pair<int, int> Bin = PhaseSpace->Bin();
      if(Bin.first != 0) {
	SignalBin = Bin.first;
	TagBin = Bin.second;
      } else {
	EventsOutsidePhaseSpace++;
	continue;
      }
    }
    if(settings.getB("Bin_truth")) {
      std::pair<int, int> Bin = PhaseSpace->TrueBin();
      if(Bin.first != 0) {
	SignalBin_true = Bin.first;
	TagBin_true = Bin.second;
      } else {
	EventsOutsidePhaseSpace++;
	continue;
      }
    }
    OutputTree->Fill();
  }
  if(settings.getB("Bin_reconstructed")) {
    std::cout << "Reconstructed events outside of phase space: " << EventsOutsidePhaseSpace << "\n";
  }
  if(settings.getB("Bin_truth")) {
    std::cout << "True events outside of phase space: " << EventsOutsidePhaseSpace_true << "\n";
  }
  OutputTree->Write();
  OutputFile.Close();
  std::cout << "Binning complete\n";
  return 0;
}
