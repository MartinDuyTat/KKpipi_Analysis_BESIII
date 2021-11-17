// Martin Duy Tat 28th April 2021
/**
 * DetermineDoubleTagYields is an application that calculates the double tag yields of a sample of double tagged events and saves them to a text file
 * @param 1 Filename of ROOT file with double tagged events
 * @param 2 Name of TTree with double tagged events
 * @param 3 Number of bins in the binning scheme
 * @param 4 Filename of text file where yields are saved
 * @param 5 "Flavour" or "CP"
 * @param 6 (optional) Cuts to apply to TTree before counting events
 */

#include<iostream>
#include<fstream>
#include<cstdlib>
#include"TChain.h"
#include"DoubleTagYield.h"

int main(int argc, char *argv[]) {
  if(argc != 6 && argc != 7) {
    std::cout << "Need 5 or 6 input arguments\n";
    return 0;
  }
  std::string TagType(argv[5]);
  if(TagType != "Flavour" && TagType != "CP") {
    std::cout << TagType << " is an unknown tag type\n";
    return 0;
  }
  TChain Chain(argv[2]);
  Chain.Add(argv[1]);
  int NBins = std::atoi(argv[3]);
  TCut Cuts = argc == 5 ? TCut() : TCut(argv[6]);
  DoubleTagYield doubleTagYield(NBins, Cuts);
  doubleTagYield.CalculateBinnedRawYields(&Chain);
  std::ofstream Outfile(argv[4]);
  for(int i = 1; i <= NBins; i++) {
    double BinYieldPlus = doubleTagYield.GetBinYield(i);
    double BinYieldMinus = doubleTagYield.GetBinYield(-i);
    if(TagType == "CP") {
      Outfile << "Bin " << i << ": " << BinYieldPlus + BinYieldMinus << "\n";
    } else if(TagType == "Flavour") {
      Outfile << "Bin " << i << ": " << BinYieldPlus << "\n";
      Outfile << "Bin " << -i << ": " << BinYieldMinus << "\n";
    }
  }
  Outfile << "Total yield: " << doubleTagYield.GetTotalYield() << "\n";
  Outfile << "Events outside MBC region: " << doubleTagYield.GetEventsOutsideMBCSpace() << "\n";
  Outfile << "Events outside phase space: " << doubleTagYield.GetEventsOutsidePhaseSpace() << "\n";
  Outfile.close();
  return 0;
}
