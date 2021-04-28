// Martin Duy Tat 28th April 2021
/** 
 * PredictDoubleTaggedYields is an application for prediction the total and binned yields of \f$KK\pi\pi\f$, tagged by a \f$CP = \pm\f$ or flavour tag
 * @param 1 Filename of file with \f$D\f$ hadronic parameters
 * @param 2 Double tag yield of flavour tags
 * @param 3 Single tag yield of flavour tags (not required for flavour tag yield)
 * @param 4 Single tag yield of \f$CP\f$ tag (not required for flavour tag yield)
 * @param 5 \f$CP\f$ of tag mode (not required for flavour tag yield)
 */

#include<iostream>
#include<cstdlib>
#include"PredictNumberEvents.h"
#include"DDecayParameters.h"

int main(int argc, char *argv[]) {
  if(argc != 3 && argc != 6) {
    std::cout << "Need 2 or 5 input arguments\n";
    return 0;
  }
  std::cout << "Calculation of event yields for double tagged events\n";
  double DTFlavourYield = std::atof(argv[2]);
  double STFlavourYield = argc == 3 ? 0.0 : std::atof(argv[3]);
  double STYield = argc == 3 ? 0.0 : std::atof(argv[4]);
  double CP = argc == 3 ? 0 : std::atoi(argv[5]);
  PredictNumberEvents predictNumberEvents(DDecayParameters(std::string(argv[1])), DTFlavourYield, STFlavourYield, STYield, CP);
  int NBins = predictNumberEvents.GetNBins();
  double TotalYield = 0.0;
  std::cout << "Binned event yields:\n";
  for(int i = 1; i <= NBins; i++) {
    if(i == 0) {
      continue;
    }
    double BinYieldPlus = predictNumberEvents(i);
    double BinYieldMinus = predictNumberEvents(-i);
    TotalYield += BinYieldPlus + BinYieldMinus;
    if(CP == 0) {
      std::cout << "Bin " << i << ": " << BinYieldPlus << "\n";
      std::cout << "Bin " << -i << ": " << BinYieldMinus << "\n";
    } else {
      std::cout << "Bin " << i << ": " << BinYieldPlus + BinYieldMinus << "\n";
    }
  }
  std::cout << "Total integrated yield: " << TotalYield << "\n";
  return 0;
}
