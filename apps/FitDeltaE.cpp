// Martin Duy Tat 31st March 2021
/**
 * FitDeltaE performs a fit of \f$\Delta E\f$ on an input ROOT file with a TTree
 */

#include<iostream>
#include"TFile.h"
#include"TChain.h"
#include"DeltaEFit.h"
#include"Utilities.h"
#include"Settings.h"

int main(int argc, char *argv[]) {
  std::cout << "Fit of Delta E\n";
  std::cout << "Loading settings...\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "Settings ready\n";
  std::cout << "Loading ROOT file\n";
  TChain Chain(settings.get("TreeName").c_str());
  Utilities::LoadChain(&Chain, settings.get("InputFiles"));
  DeltaEFit deltaEFit(&Chain, settings);
  deltaEFit.FitDeltaE();
  //deltaEFit.SaveParameters();
  return 0;
}
