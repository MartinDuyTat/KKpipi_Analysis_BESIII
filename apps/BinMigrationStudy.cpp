// Martin Duy Tat 3rd December 2021
/**
 * BinMigrationStudy is an application that takes in a set of D daughter momenta and smears them before determining the amplitude ratio \f$r_D\f$ and phase difference \f$\delta_D\f$ using the LHCb amplitude model
 * The momenta are space-separated numbers read line by line from a file, and for each line a new TTree is created
 */

#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<complex>
#include"TTree.h"
#include"TFile.h"
#include"TH1D.h"
#include"TRandom.h"
#include"Settings.h"
#include"Utilities.h"
#include"Amplitude.h"

int main(int argc, char *argv[]) {
  std::cout << "Bin migration study\n";
  Settings settings = Utilities::parse_args(argc, argv);
  // Vectors defining the order of the particles
  std::vector<std::string> Particles{"Kplus", "Kminus", "piplus", "piminus"};
  std::vector<std::string> Components{"PX", "PY", "PZ", "PE"};
  std::cout << "Getting the resolution from the resolution histograms...\n";
  std::vector<double> Resolutions;
  std::string ResolutionFilename = settings.get("ResolutionFilename");
  TFile ResolutionFile(ResolutionFilename.c_str(), "READ");
  for(const auto &Particle : Particles) {
    for(const auto &Component : Components) {
      TH1D *Histogram = nullptr;
      ResolutionFile.GetObject((Particle + Component + "_h").c_str(), Histogram);
      Resolutions.push_back(Histogram->GetStdDev());
    }
  }
  ResolutionFile.Close();
  std::cout << "Resolutions loaded\n";
  std::cout << "Preparing input and output files...\n";
  std::string OutputFilename = settings.get("OutputFilename");
  TFile OutputFile(OutputFilename.c_str(), "RECREATE");
  // Open file with initial D daughter momenta
  std::string DaughterMomentaFilename = settings.get("DaughterMomentaFilename");
  std::ifstream DaughterMomentaFile(DaughterMomentaFilename);
  std::cout << "Files open and ready\n";
  std::string Line;
  // Prepare phase space
  Amplitude amplitude;
  int TreeNumber = 0;
  while(std::getline(DaughterMomentaFile, Line)) {
    std::cout << "Smearing point number " << TreeNumber << "\n";
    // Parse momentum components
    std::stringstream ss(Line);
    std::vector<double> Momentum(16);
    for(int i = 0; i < 16; i++) {
      ss >> Momentum[i];
    }
    // Create TTree containing the smeared components
    TTree Tree(("BinMigrationTree" + std::to_string(TreeNumber)).c_str(), "");
    double rD, deltaD;
    Tree.Branch("rD", &rD);
    Tree.Branch("deltaD", &deltaD);
    // First entry is without smearing
    std::complex<double> Amp_D0 = amplitude(Momentum, +1);
    std::complex<double> Amp_D0bar = amplitude(Momentum, -1);
    rD = std::abs(Amp_D0/Amp_D0bar);
    deltaD = std::arg(Amp_D0/Amp_D0bar);
    Tree.Fill();
    for(int i = 0; i < settings.getI("NumberSmearings"); i++) {
      if(i%(settings.getI("NumberSmearings")/10) == 0) {
	std::cout << "Smearing " << i << "\n";
      }
      std::vector<double> SmearedMomenta = Momentum;
      for(int j = 0; j < 16; j++) {
	SmearedMomenta[j] += gRandom->Gaus(0.0, Resolutions[j]);
      }
      Amp_D0 = amplitude(SmearedMomenta, +1);
      amplitude(SmearedMomenta, -1);
      rD = std::abs(Amp_D0/Amp_D0bar);
      deltaD = std::arg(Amp_D0/Amp_D0bar);
      Tree.Fill();
    }
    TreeNumber++;
    Tree.Write();
  }
  std::cout << "Bin migration study completed\n";
  return 0;
}
