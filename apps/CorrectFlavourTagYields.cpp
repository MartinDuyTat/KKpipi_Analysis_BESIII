// Martin Duy Tat 1st December 2021
/**
 * Correct flavour tag yields is an application that performs efficiency and DCS corrections to the \f$K_i\f$
 */

#include<iostream>
#include<fstream>
#include"TMatrixT.h"
#include"TFile.h"
#include"TMath.h"
#include"Settings.h"
#include"Utilities.h"
#include"Category.h"

int main(int argc, char *argv[]) {
  std::cout << "Making efficiency and DCS corrections to flavour tag yields\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::string Mode = settings.get("Mode");
  int NumberBins = settings["BinningScheme"].getI("NumberBins");
  Category category(settings);
  std::cout << "Loading double tag yields...\n";
  TMatrixT<double> Yields(2*NumberBins, 1);
  TMatrixT<double> YieldErrors(2*NumberBins, 1);
  int i = 0;
  for(const auto &cat : category.GetCategories()) {
    Yields(i, 0) = settings["DoubleTagYields"].getD(cat + "_SignalYield");
    YieldErrors(i, 0) = settings["DoubleTagYields"].getD(cat + "_SignalYield_err");
    i++;
  }
  std::cout << "Yields ready\n";
  std::cout << "Loading efficiency matrix...\n";
  TFile EffMatrixFile(settings.get("EfficiencyMatrixFile").c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  EffMatrixFile.GetObject("EffMatrix", EffMatrix);
  EffMatrix->Invert();
  TMatrixT<double> EffCorrectedYields = Yields;
  EffCorrectedYields = *EffMatrix*EffCorrectedYields;
  std::cout << "Yields efficiency corrected\n";
  std::cout << "Making DCS corrections...\n";
  std::ofstream Outfile(settings.get("Ki_Filename"));
  i = 0;
  double Sum = 0.0;
  // Ignore uncertainties from efficiency and DCS corrections for now, just scale them
  for(int Bin = -NumberBins; Bin <= NumberBins; Bin++) {
    if(Bin == 0) {
      continue;
    }
    std::string DCS_Name = Mode + "_DCS_Correction_";
    DCS_Name += Bin > 0 ? "P" : "M";
    DCS_Name += std::to_string(TMath::Abs(Bin));
    EffCorrectedYields(i, 0) *= settings["DCS_Corrections"].getD(DCS_Name);
    std::string VariableName = Mode + "_Ki_Bin" + (Bin > 0 ? "P" : "M" + std::to_string(TMath::Abs(Bin)));
    Outfile << VariableName << " " << EffCorrectedYields(i, 0) << "\n";
    YieldErrors(i, 0) = EffCorrectedYields(i, 0)*YieldErrors(i, 0)/Yields(i, 0);
    Outfile << VariableName << "_err " << YieldErrors(i, 0) << "\n";
    Sum += EffCorrectedYields(i, 0);
    i++;
  }
  Outfile << "\n";
  i = 0;
  for(int Bin = -NumberBins; Bin <= NumberBins; Bin++) {
    if(Bin == 0) {
      continue;
    }
    std::string VariableName = Mode + "_Ki_Norm_Bin" + (Bin > 0 ? "P" : "M" + std::to_string(TMath::Abs(Bin)));
    Outfile << VariableName << " " << EffCorrectedYields(i, 0)/Sum << "\n";
    Outfile << VariableName << "_err " << YieldErrors(i, 0)/Sum << "\n";
    i++;
  }
  Outfile.close();
  std::cout << "Yields are now efficiency and DCS corrected!\n";
  return 0;
}
