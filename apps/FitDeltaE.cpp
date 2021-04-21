// Martin Duy Tat 31st March 2021
/**
 * FitDeltaE performs a fit of \f$\Delta E\f$ on an input ROOT file with a TTree
 * @param 1 Filename of text files with a list of ROOT files with TTree objects
 * @param 2 Name of TTree
 * @param 3 Name of text file where parameters are saved
 * @param 4 Tag mode label in plots
 * @param 5 Filename for plots
 * @param 6 Input "Unbinned" to do an unbinned fit after the binned fit, "Binned" for a binned fit only and "NoFit" if no fit is performed at all
 * @param 7 Filename of file with alternative parameters for the fit (optional)
 */

#include<iostream>
#include"TFile.h"
#include"TChain.h"
#include"DeltaEFit.h"
#include"Utilities.h"

int main(int argc, char *argv[]) {
  if(argc != 7 && argc != 8) {
    std::cout << "Need 6 or 7 input arguments\n";
    return 0;
  }
  std::cout << "Fit of Delta E\n";
  std::cout << "Loading ROOT file\n";
  TChain Chain(argv[2]);
  Utilities::LoadChain(&Chain, std::string(argv[1]));
  DeltaEFit deltaEFit(&Chain);
  bool DoBinnedFit = false;
  bool DoUnbinnedFit = false;
  if(std::string(argv[6]) == "Binned") {
    DoBinnedFit = true;
  } else if(std::string(argv[6]) == "Unbinned") {
    DoBinnedFit = true;
    DoUnbinnedFit = true;
  }
  if(argc == 8) {
    deltaEFit.SetParameters(std::string(argv[7]));
  }
  deltaEFit.FitDeltaE(std::string(argv[5]), std::string(argv[4]), DoBinnedFit, DoUnbinnedFit);
  deltaEFit.SaveParameters(std::string(argv[3]));
  return 0;
}
