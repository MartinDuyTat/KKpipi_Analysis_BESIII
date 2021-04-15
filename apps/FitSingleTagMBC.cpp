// Martin Duy Tat 4th April 2021
/**
 * FitSingleTagMBC is an application that determines the single tag yield using a fit of the beam constrained mass \f$m_\text{BC}\f$ distribution
 * The signal shape is taken from an exclusive signal MC sample, convolved with a Gaussian, and the background is modelled with an Argus PDF
 * @param 1 Filename of text file with filename of data samples
 * @param 2 Filename of ROOT file with signal MC sample
 * @param 3 Name of TTree
 * @param 4 Tag mode label for plots
 * @param 5 Filename to save plot to
 * @param 6 Filename to save fit parameters to
 * @param 7 Filename with peaking background components (optional)
 */

#include<iostream>
#include<string>
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"SingleTagYield.h"
#include"Utilities.h"

int main(int argc, char *argv[]) {
  if(argc != 7 && argc != 8) {
    std::cout << "Need 6 or 7 input arguments\n";
    return 0;
  }
  std::cout << "Single tag yield fit\n";
  std::cout << "Loading ROOT files...\n";
  TChain Chain(argv[3]);
  Utilities::LoadChain(&Chain, std::string(argv[1]));
  TFile MCFile(argv[2], "READ");
  TTree *MCTree = nullptr;
  MCFile.GetObject(argv[3], MCTree);
  std::cout << "Trees ready\n";
  std::cout << "Fitting yield...\n";
  SingleTagYield singleTagYield(&Chain, MCTree);
  if(argc == 8) {
    singleTagYield.AddPeakingComponent(std::string(argv[7]));
  }
  singleTagYield.FitYield(std::string(argv[4]), std::string(argv[5]));
  singleTagYield.SaveFitParameters(std::string(argv[6]));
  std::cout << "Fit completed, yield is " << singleTagYield.GetYield() << "\n";
  return 0;
}
