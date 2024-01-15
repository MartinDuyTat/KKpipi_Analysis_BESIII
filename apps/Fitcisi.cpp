// Martin Duy Tat 13th October 2022
/**
 * Fitcisi is an application that loads double tag yields, single tag yields, efficiencies and fits the strong phases \f$c_i\f$ and \f$s_i\f$
 */

#include<iostream>
#include<string>
#include<vector>
#include"Settings.h"
#include"Utilities.h"
#include"cisiFitter.h"

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "ci and si fit\n";
  std::vector<std::string> TagModes =
    Utilities::ConvertStringToVector(settings.get("TagModes"));
  cisiFitter Fitter(settings);
  const std::string RunMode(settings.get("RunMode"));
  if(RunMode == "SingleFit") {
    std::cout << "Fitting ci and si...\n";
    Fitter.Minimise();
    std::cout << "ci and si have been measured!\n";
  } else if(RunMode == "SingleToy") {
    std::cout << "Running ci and si toys...\n";
    int ToyNumber = settings.getI("ToyNumber");
    if(ToyNumber < 0) {
      return 0;
    }
    Fitter.RunToy(ToyNumber);
    std::cout << "Toy studies of ci and si are ready!\n";
  } else if(RunMode == "DumpGeneratorYields") {
    std::cout << "Dumping generator yields for toy studies...\n";
    Fitter.SavePredictedYields();
    std::cout << "Ready to do toys!\n";
  } else if(RunMode == "DumpGeneratorYieldsForFeldmanCousins") {
    std::cout << "Dumping generator yields for a Feldman Cousins scan...\n";
    std::size_t Parameter = settings.getI("FeldmanCousins_Parameter");
    Fitter.GenerateFeldmanCousinsYields(Parameter);
    std::cout << "Ready to do Feldman Cousins scan!\n";
  } else if(RunMode == "FitFeldmanCousinsToy") {
    std::string ToyName = settings.get("FeldmanCousins_ToyName");
    if(ToyName == "None") {
      return 0;
    }
    std::cout << "Fitting Feldman Cousins toy " << ToyName;
    int ToyNumber = settings.getI("ToyNumber");
    if(ToyNumber < 0) {
      return 0;
    }
    int Parameter = std::stoi(ToyName.substr(13, 1));
    Fitter.FitFeldmanCousinsToy(ToyName, ToyNumber, Parameter);
    std::cout << "Feldman Cousins Delta chi2 saved!\n";
  } else if(RunMode == "FeldmanCousinsDataScan") {
    std::cout << "Doing a Feldman Cousins scan over data...\n";
    Fitter.FeldmanCousinsDataScan();
    std::cout << "Feldman Cousins scan complete!\n";
  } else if(RunMode == "SingleFitWithGamma") {
    std::cout << "Fitting ci and si simultaneously with gamma...\n";
    Fitter.MinimiseWithGamma();
    std::cout << "Model-independent measurement of gamma complete!\n";
  } else if(RunMode == "RunFeldmanCousinsWithGamma") {
    const double gamma = settings.getD("FeldmanCousins_gamma_value");
    std::cout << "Doing a Feldman Cousins gamma scan with gamma = ";
    std::cout << gamma << "\n";
    if(gamma < 0.0 || gamma > 180.0) {
      std::cout << "gamma must be in the range [0, 180] degrees!\n";
      return 0;
    }
    Fitter.ScanGammaFeldmanCousins(gamma);
    std::cout << "Feldman Cousins scan with gamma complete!\n";
  } else if(RunMode == "ManyGammaToys") {
    std::cout << "Running toy studies with combined fit of gamma, ci and si...\n";
    Fitter.RunGammaToys();
    std::cout << "Toy studies of gamma, ci and si are ready!\n";
  }
  return 0;
}
