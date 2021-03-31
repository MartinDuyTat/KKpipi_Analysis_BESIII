// Martin Duy Tat 31st March 2021
/**
 * FitDeltaE performs a fit of \f$\Delta E\f$ on an input ROOT file with a TTree
 * @param 1 Filename of ROOT file with TTree
 * @param 2 Name of TTree
 * @param 3 Name of text file where parameters are saved
 * @param 4 Tag mode label in plots
 * @param 5 Filename for plots
 * @param 6 Input "yes" to do an unbinned fit after the binned fit
 */

#include<iostream>
#include"TFile.h"
#include"TTree.h"
#include"DeltaEFit.h"

int main(int argc, char *argv[]) {
  if(argc != 6 && argc != 7) {
    std::cout << "Need 5 or 6 input arguments\n";
    return 0;
  }
  std::cout << "Fit of Delta E\n";
  std::cout << "Loading ROOT file\n";
  TFile Infile(argv[1], "READ");
  TTree *Tree = nullptr;
  Infile.GetObject(argv[2], Tree);
  std::cout << "TTree ready to be analyzed\n";
  DeltaEFit deltaEFit(Tree);
  bool DoUnbinnedFit = argc == 7 && std::string(argv[6]) == "yes";
  deltaEFit.FitDeltaE(std::string(argv[5]), std::string(argv[4]), DoUnbinnedFit);
  deltaEFit.SaveParameters(std::string(argv[3]));
  Infile.Close();
  return 0;
}
