// Martin Duy Tat 3rd December 2021
/**
 * MakeResolutionHistograms will loop over an MC sample and create a histogram for each \f$D\f$ daughter in the \f$KK\pi\pi\f$ decay
 */

#include<iostream>
#include<string>
#include<map>
#include<vector>
#include"TChain.h"
#include"TH1D.h"
#include"TFile.h"
#include"Utilities.h"
#include"Settings.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"

int main(int argc, char *argv[]) {
  std::cout << "Making D->KKpipi daughter resolution histograms\n";
  Settings settings = Utilities::parse_args(argc, argv);
  std::string TreeName = settings.get("TreeName");
  std::cout << "Loading double tag events and setting up output TTree...\n";
  TChain InputChain(TreeName.c_str());
  InputChain.Add(settings.get("InputFilename").c_str());
  std::unique_ptr<KKpipi_PhaseSpace> PhaseSpace = Utilities::GetPhaseSpaceBinning(settings, &InputChain);
  std::cout << "TTree and phase space ready\n";
  std::string HistogramsFilename = settings.get("HistogramsFilename");
  TFile OutputFile(HistogramsFilename.c_str(), "RECREATE");
  std::vector<std::string> Particles{"Kplus", "Kminus", "piplus", "piminus"};
  std::vector<std::string> Components{"PX", "PY", "PZ", "PE"};
  std::map<std::string, TH1D> Histograms;
  std::cout << "Filling histograms...\n";
  for(const auto &Particle : Particles) {
    for(const auto &Component : Components) {
      Histograms.insert({Particle + Component, TH1D((Particle + Component + "_h").c_str(), "", 100, 0.0, 0.0)});
    }
  }
  for(int i = 0; i < InputChain.GetEntries(); i++) {
    InputChain.GetEntry(i);
    PhaseSpace->TrueBin();
    auto Resolution = PhaseSpace->GetMomentumResolution();
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
	Histograms.at(Particles[j] + Components[k]).Fill(Resolution[4*j + k]);
      }
    }
  }
  OutputFile.Write();
  OutputFile.Close();
  std::cout << "Resolution histograms saved\n";
  return 0;
}
