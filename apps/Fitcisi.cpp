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
  std::cout << "Fitting ci and si...\n";
  Fitter.Minimise();
  std::cout << "ci and si have been measured!\n";
  return 0;
}
