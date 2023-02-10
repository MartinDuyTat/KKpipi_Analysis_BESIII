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
  } else if(RunMode == "ManyToys") {
    std::cout << "Running ci and si toys...\n";
    Fitter.RunToys();
    std::cout << "Toy studies of ci and si are ready!\n";
  }
  return 0;
}
