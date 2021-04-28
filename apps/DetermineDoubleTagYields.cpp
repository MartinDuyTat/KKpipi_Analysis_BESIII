// Martin Duy Tat 28th April 2021
/**
 * DetermineDoubleTagYields is an application that calculates the double tag yields of a sample of double tagged events and saves them to a text file
 * @param 1 Filename of ROOT file with double tagged events
 * @param 2 Name of TTree with double tagged events
 * @param 3 Number of bins in the binning scheme
 * @param 4 Filename of text file where yields are saved
 * @param 5 (optional) Cuts to apply to TTree before counting events
 */

#include<iostream>
#include<fstream>
#include<cstdlib>
#include"TChain.h"
#include"DoubleTagYield.h"

int main(int argc, char *argv[]) {
  if(argc != 5 && argc != 6) {
    std::cout << "Need 4 or 5 input arguments\n";
    return 0;
  }
  TChain Chain(argv[2]);
  Chain.Add(argv[1]);
  int NBins = std::atoi(argv[3]);
  TCut Cuts = argc == 5 ? TCut() : TCut(argv[5]);
  DoubleTagYield doubleTagYield(NBins, Cuts);
  doubleTagYield.CalculateBinnedRawYields(&Chain);
  std::ofstream Outfile(argv[4]);
  double TotalYield = 0.0;
  for(int i = 1; i <= NBins; i++) {
    double BinYield = doubleTagYield.GetBinYield(i);
    TotalYield += BinYield > 0.0 ? BinYield : 0.0;
    Outfile << "Bin " << i << ": " << BinYield << "\n";
  }
  Outfile << "Total yield: " << TotalYield << "\n";
  Outfile << "Events outside MBC region: " << doubleTagYield.GetEventsOutsideMBCSpace() << "\n";
  Outfile << "Events outside phase space: " << doubleTagYield.GetEventsOutsidePhaseSpace() << "\n";
  Outfile.close();
  return 0;
}
