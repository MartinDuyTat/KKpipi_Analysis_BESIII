// Martin Duy Tat 28th April 2021
/**
 * FitFPlus is an application that loads double tag yields, normalized by their single tag yields, and performs a fit to determine the CP even fraction \f$F_+\f$
 */

#include<iostream>
#include<string>
#include<vector>
#include"Settings.h"
#include"Utilities.h"
#include"FPlusFitter.h"

int main(int argc, char *argv[]) {
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "F+ fit\n";
  std::vector<std::string> TagModes = Utilities::ConvertStringToVector(settings.get("TagModes"));
  FPlusFitter Fitter(settings);
  for(const auto &TagMode : TagModes) {
    Fitter.AddTag(TagMode);
  }
  std::cout << "Fitting F+...\n";
  Fitter.InitializeAndFit();
  std::cout << "F+ has been measured!\n";
  return 0;
}
