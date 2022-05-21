// Martin Duy Tat 26th November 2021
/**
 * BinDoubleTags is an application that determines which phase space bin a double tag belong to
 * Both the signal and tag side are analyzed and a combined phase space bin is saved to the ROOT file
 */

#include<iostream>
#include<memory>
#include<string>
#include<utility>
#include<stdexcept>
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
  std::vector<std::string> DalitzVariables{"s01", "s03", "s12", "s23", "s012"};
  std::map<std::string, double> DalitzCoordinates, RecDalitzCoordinates;
  TFile OutputFile(OutputFilename.c_str(), "RECREATE");
  TTree *OutputTree = InputChain.CloneTree(0);
  if(settings.getB("Bin_reconstructed")) {
    OutputTree->Branch(SignalBin_Name.c_str(), &SignalBin);
    OutputTree->Branch(TagBin_Name.c_str(), &TagBin);
    for(const auto &DalitzVariable : DalitzVariables) {
      RecDalitzCoordinates.insert({DalitzVariable, 0.0});
      OutputTree->Branch(("Rec" + DalitzVariable).c_str(), &RecDalitzCoordinates[DalitzVariable]);
    }
  }
  if(settings.getB("Bin_truth")) {
    OutputTree->Branch((SignalBin_Name + "_true").c_str(), &SignalBin_true);
    OutputTree->Branch((TagBin_Name + "_true").c_str(), &TagBin_true);
    for(const auto &DalitzVariable : DalitzVariables) {
      DalitzCoordinates.insert({DalitzVariable, 0.0});
      OutputTree->Branch(DalitzVariable.c_str(), &DalitzCoordinates[DalitzVariable]);
    }
  }
  int EventsOutsidePhaseSpace = 0;
  int EventsOutsidePhaseSpace_true = 0;
  int NumberExceptions = 0;
  std::cout << "Ready to bin phase space\n";
  int Entries = InputChain.GetEntries();
  for(int i = 0; i < Entries; i++) {
    InputChain.GetEntry(i);
    if(settings.getB("Bin_reconstructed")) {
      std::pair<int, int> Bin = PhaseSpace->Bin();
      SignalBin = Bin.first;
      TagBin = Bin.second;
      if(Bin.first == 0) {
	EventsOutsidePhaseSpace++;
	if(!(settings.contains("IncludeEventsOutsidePhaseSpace") && 
	     settings.getB("IncludeEventsOutsidePhaseSpace"))) {
	  continue;
	}
      }
      auto RecoDalitzCoordinates = PhaseSpace->GetRecDalitzCoordinates();
      for(const auto &DalitzVariable : DalitzVariables) {
	RecDalitzCoordinates[DalitzVariable] = RecoDalitzCoordinates[DalitzVariable];
      }
    }
    if(settings.getB("Bin_truth")) {
      std::pair<int, int> Bin;
      try {
	Bin = PhaseSpace->TrueBin();
      } catch(const std::logic_error &e) {
	NumberExceptions++;
	  continue;
      }
      SignalBin_true = Bin.first;
      TagBin_true = Bin.second;
      if(Bin.first == 0) {
	EventsOutsidePhaseSpace_true++;
	if(!(settings.contains("IncludeEventsOutsidePhaseSpace") && 
	     settings.getB("IncludeEventsOutsidePhaseSpace"))) {
	  continue;
	}
      }
      auto TrueDalitzCoordinates = PhaseSpace->GetDalitzCoordinates();
      for(const auto &DalitzVariable : DalitzVariables) {
	DalitzCoordinates[DalitzVariable] = TrueDalitzCoordinates[DalitzVariable];
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
  std::cout << "Number of events caught and elegantly skipped: " << NumberExceptions << "\n";
  OutputTree->Write();
  OutputFile.Close();
  std::cout << "Binning complete\n";
  return 0;
}
