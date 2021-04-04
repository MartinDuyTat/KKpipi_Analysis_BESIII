// Martin Duy Tat 4th April 2021
/**
 * FitSingleTagMBC is an application that determines the single tag yield using a fit of the beam constrained mass \f$m_\text{BC}\f$ distribution
 * The signal shape is taken from an exclusive signal MC sample, convolved with a Gaussian, and the background is modelled with an Argus PDF
 * @param 1 Filename of ROOT file with data sample
 * @param 2 Filename of ROOT file with signal MC sample
 * @param 3 Name of TTree
 * @param 4 Tag mode label for plots
 * @param 5 Filename to save plot to
 * @param 6 Filename to save fit parameters to
 */

#include<iostream>
#include<string>
#include"TFile.h"
#include"TTree.h"
#include"SingleTagYield.h"

int main(int argc, char *argv[]) {
  if(argc != 7) {
    std::cout << "Need 6 input arguments\n";
    return 0;
  }
  std::cout << "Single tag yield fit\n";
  std::cout << "Loading ROOT files...\n";
  TFile DataFile(argv[1], "READ");
  TTree *DataTree = nullptr;
  DataFile.GetObject(argv[3], DataTree);
  TFile MCFile(argv[2], "READ");
  TTree *MCTree = nullptr;
  MCFile.GetObject(argv[3], MCTree);
  std::cout << "Trees ready\n";
  std::cout << "Fitting yield...\n";
  SingleTagYield singleTagYield(DataTree, MCTree);
  singleTagYield.FitYield(std::string(argv[4]), std::string(argv[5]));
  singleTagYield.SaveFitParameters(std::string(argv[6]));
  std::cout << "Fit completed, yield is " << singleTagYield.GetYield() << "\n";
  return 0;
}
