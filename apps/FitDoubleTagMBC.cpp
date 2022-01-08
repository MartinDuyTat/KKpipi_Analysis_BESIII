// Martin Duy Tat 28th April 2021
/**
 * DetermineDoubleTagYields is an application that calculates the double tag yields of a sample of double tagged events and saves them to a text file
 */

#include<iostream>
#include<string>
#include"TChain.h"
#include"DoubleTagYield.h"
#include"Settings.h"
#include"Utilities.h"

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "Double tag yield fit\n";
  std::cout << "Loading ROOT files...\n";
  std::string TreeName = settings.get("TreeName");
  TChain Chain(TreeName.c_str());
  std::string Mode = settings.get("Mode");
  std::string Filename = Utilities::ReplaceString(settings["BinnedDataSets"].get("BinnedDataSet").c_str(), "TAG", Mode);
  Chain.Add(Filename.c_str());
  DoubleTagYield doubleTagYield(settings, &Chain);
  doubleTagYield.DoFit();
  return 0;
}
