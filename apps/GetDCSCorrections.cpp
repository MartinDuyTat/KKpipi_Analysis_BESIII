// Martin Duy Tat 1st December 2021
/**
 * GetDCSCorrections is an application that calculates the DCS corrections for the flavour tags
 */

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include"TMath.h"
#include"Utilities.h"
#include"uncertainties/ureal.hpp"
#include"uncertainties/math.hpp"
#include"HadronicParameters/cisi.h"
#include"HadronicParameters/Ki.h"
#include"HadronicParameters/DCS_Parameters.h"

int main(int argc, char *argv[]) {
  std::cout << "Calculating DCS corrections\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::string Mode = settings.get("Mode");
  std::cout << "Initializing hadronic parameters...\n";
  cisi cisi_m(settings["BinningScheme"]);
  Ki Ki_m(settings["BinningScheme"]);
  DCS_Parameters DCS(settings["DCS_Parameters_" + Mode]);
  std::cout << "Hadronic parameters ready\n";
  std::vector<uncertainties::udouble> DCS_parameters = DCS.GetDCSParameters();
  std::ofstream Outfile(settings.get("DCS_Parameters_Filename"));
  int NumberBins = settings["BinningScheme"].getI("NumberBins");
  std::cout << "Calculating DCS corrections...\n";
  for(int Bin = -NumberBins; Bin <= NumberBins; Bin++) {
    if(Bin == 0) {
      continue;
    }
    uncertainties::udouble Correction = Ki_m.Get_Ki(Bin)*Ki_m.Get_Ki(Bin)/(Ki_m.Get_Ki(Bin)*Ki_m.Get_Ki(Bin) + DCS_parameters[0]*DCS_parameters[0]*Ki_m.Get_Ki(-Bin)*Ki_m.Get_Ki(-Bin) - 2*DCS_parameters[0]*DCS_parameters[1]*Ki_m.Get_Ki(Bin)*Ki_m.Get_Ki(-Bin)*(cisi_m.Get_ci(Bin)*uncertainties::cos(TMath::Pi()*DCS_parameters[2]/180.0) - cisi_m.Get_si(Bin)*uncertainties::sin(TMath::Pi()*DCS_parameters[2]/180.0)));
    std::string BinSign = Bin > 0 ? "P" : "M";
    Outfile << Mode << "_DCS_Correction_" << BinSign << TMath::Abs(Bin) << " " << uncertainties::nom(Correction) << "\n";
    Outfile << Mode << "_DCS_Correction_" << BinSign << TMath::Abs(Bin) << "_err " << uncertainties::sdev(Correction) << "\n";
  }
  std::cout << "DCS corrections calculated!\n";
  return 0;
}
