// Martin Duy Tat 26th November 2021
/**
 * GetSingleTagEfficiencies is an application that calculates all the single tag efficiencies using signal MC samples
 */

#include<iostream>
#include<fstream>
#include<string>
#include"TChain.h"
#include"TMath.h"
#include"Utilities.h"
#include"Settings.h"

int main(int argc, char *argv[]) {
  std::cout << "Calculating single tag efficiencies from signal MC\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::ofstream Outfile(settings.get("SingleTagEfficienciesFilename"));
  Outfile << "* Single tag efficiencies\n\n";
  double SignalMC_SampleSize = static_cast<double>(settings["SignalMCSampleSize"].getI("SingleTag"));
  std::string ModeList = settings.get("ListModes");
  std::vector<std::string> Modes = Utilities::ConvertStringToVector(ModeList);
  for(const auto &Mode : Modes) {
    std::cout << "Analyzing " << Mode << " single tag efficiency...\n";
    std::string TreeName = settings.get("TreeName");
    TreeName = Utilities::ReplaceString(TreeName, "TAG", Mode);
    TChain Chain(TreeName.c_str());
    std::string Filename = settings["Datasets_WithDeltaECuts"].get("SignalMC_ST");
    Filename = Utilities::ReplaceString(Filename, "TAG", Mode);
    Chain.Add(Filename.c_str());
    double Entries = static_cast<double>(Chain.GetEntries("MBC > 1.86 && MBC < 1.87"));
    double Efficiency = Entries/SignalMC_SampleSize;
    Outfile << Mode << "_SingleTagEfficiency     " << Efficiency << "\n";
    Outfile << Mode << "_SingleTagEfficiency_err " << TMath::Sqrt(Efficiency*(1 - Efficiency)/SignalMC_SampleSize) << "\n";
  }
  Outfile.close();
  std::cout << "Single tag efficiencies saved\n";
  return 0;
}
