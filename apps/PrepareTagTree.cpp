// Martin Duy Tat 31st March 2021
/**
 * PrepareTagTree is an application that takes in BESIII events selected in BOSS and applies initial cuts to a TTree, which is saved to a separate file
 */

#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<sstream>
#include"TFile.h"
#include"TChain.h"
#include"TDirectory.h"
#include"TTree.h"
#include"Utilities.h"
#include"ApplyCuts.h"
#include"Settings.h"

std::vector<std::string> ParseDatasets(std::string DatasetsString);

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::string TagType = settings.get("TagType");
  std::string SignalMode("");
  if(TagType == "DT") {
    SignalMode = settings.get("SignalMode");
  }
  std::string TagMode = settings.get("TagMode");
  std::string RecSignalMode(""), RecTagMode("");
  if(settings.contains("ReconstructedSignalMode")) {
    RecSignalMode = "_to_" + settings.get("ReconstructedSignalMode");
  }
  if(settings.contains("ReconstructedTagMode")) {
    RecTagMode = "_to_" + settings.get("ReconstructedTagMode");
  }
  bool IncludeDeltaECuts = settings.getB("Include_DeltaE_Cuts");
  bool TruthMatch = settings.getB("TruthMatch");
  std::string TreeName = settings.get("TreeName");
  std::vector<std::string> Datasets = Utilities::ConvertStringToVector(settings.get("Datasets_to_include"));
  std::cout << "Sample prepration of " << TagType << " ";
  if(TagType == "DT") {
    std::cout << SignalMode << "vs ";
  }
  std::cout << TagMode << "\n";
  for(const auto &Dataset : Datasets) {
    int DataSetType = settings["DataTypes"].getI(Dataset);
    std::cout << "Processing " << Dataset << " sample, dataset type " << DataSetType << "\n";
    std::string InputFilename;
    if(DataSetType != 10) {
      InputFilename = settings["Datasets"].get(Dataset);
    } else {
      InputFilename = settings["Datasets"].get(Dataset + "_" + TagType);
      InputFilename = Utilities::ReplaceString(InputFilename, "SIGNAL", SignalMode);
      if(TagMode != "KSpipipi0") {
	InputFilename = Utilities::ReplaceString(InputFilename, "TAG", TagMode);
      } else {
	InputFilename = Utilities::ReplaceString(InputFilename, "TAG", "KSomegapipipi0");
      }
    }
    std::vector<std::string> Years = InputFilename.find("YEAR") == std::string::npos ? std::vector<std::string>{""} : std::vector<std::string>{"2010", "2011"};
    for(const auto &Year : Years) {
      double LuminosityScale = settings["LuminosityScale"].getD(Dataset + Year);
      std::cout << "Luminosity scale is " << LuminosityScale << "\n";
      std::string OutputFilename = settings.get("OutputFilenamePrefix") + "_" + Dataset + Year + ".root";
      std::string DataMC = DataSetType == 0 ? "Data" : "MC";
      TCut Cuts = Utilities::LoadCuts(SignalMode + RecSignalMode, TagMode + RecTagMode, TagType, DataMC, IncludeDeltaECuts, TruthMatch);
      // This cut removes empty NTuples
      Cuts = Cuts && TCut("!(Run == 0 && Event == 0)");
      std::cout << "Cuts ready, will apply the following cuts:\n" << Cuts.GetTitle() << "\n";
      ApplyCuts applyCuts(Cuts);
      std::cout << "Cuts ready\n";
      std::cout << "Loading TChain...\n";
      TChain Chain(TreeName.c_str());
      std::string NewFilename(InputFilename);
      if(Year != "") {
	NewFilename.replace(NewFilename.find("YEAR"), 4, Year);
      }
      Chain.Add(NewFilename.c_str());
      std::cout << "Applying cuts...\n";
      TFile OutputFile(OutputFilename.c_str(), "RECREATE");
      TTree *OutputTree = applyCuts(&Chain, DataSetType, LuminosityScale);
      OutputTree->SetDirectory(&OutputFile);
      OutputTree->Write();
      OutputFile.Close();
      // I think this line prevents a seg fault for some reason
      gDirectory->Clear();
      std::cout << "Cuts applied and events saved to file " << OutputFilename << "\n";
    }
  }
  return 0;
}
