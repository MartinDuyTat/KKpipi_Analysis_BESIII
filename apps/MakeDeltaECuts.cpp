// Martin Duy Tat 2nd April 2021
/**
 * MakeDeltaECuts is an application that will take in a file with fit parameters from a \f$\Delta E\f$ fit and calculate the \f$\Delta\f$ cut boundaries
 * The cuts are saved in the DeltaECuts directory, which will later be used to prepare samples for analysis
 * @param 1 Filename of text file with fit parameters
 * @param 2 Tag mode
 * @param 3 "MC" or "Data"
 * @param 4 Type "pi0" to get an asymmetric cut for \f$\pi^0\f$ modes
 */

#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include"TMath.h"

int main(int argc, char *argv[]) {
  if(argc != 4 && argc != 5) {
    std::cout << "Need 3 or 4 input arguments\n";
    return 0;
  }
  if(std::string(argv[3]) != "MC" && std::string(argv[3]) != "Data") {
    std::cout << "Parameter 3 must be MC or Data\n";
    return 0;
  }
  std::cout << "Calculation of deltaE cuts\n";
  std::ifstream ParameterFile(argv[1]);
  double Nsig1, Nsig2, Mean1, Mean2, Sigma1, Sigma2;
  std::string line;
  if(ParameterFile.is_open()) {
    while(std::getline(ParameterFile, line)) {
      std::stringstream ss(line);
      std::string Parameter;
      ss >> Parameter;
      if(Parameter == "Nsig1") {
	ss >> Nsig1;
      } else if(Parameter == "Nsig2") {
	ss >> Nsig2;
      } else if(Parameter == "Mean1") {
	ss >> Mean1;
      } else if(Parameter == "Mean2") {
	ss >> Mean2;
      } else if(Parameter == "Sigma1") {
	ss >> Sigma1;
      } else if(Parameter == "Sigma2") {
	ss >> Sigma2;
      }
    }
  } else {
    std::cout << argv[1] << " does not exist\n";
  }
  ParameterFile.close();
  double f = Nsig1/(Nsig1 + Nsig2);
  double Mean = f*Mean1 + (1 - f)*Mean2;
  double Sigma = TMath::Sqrt(f*Sigma1*Sigma1 + (1 - f)*Sigma2*Sigma2 + f*(1 - f)*TMath::Power(Mean1 - Mean2, 2));
  double Upper = Mean + 3*Sigma;
  double Lower = Mean - 3*Sigma;
  if(argc == 5 && std::string(argv[4]) == "pi0") {
    Lower = Mean - 4*Sigma;
  }
  std::ofstream DeltaEFile(std::string(DELTAE_CUTS_DIR) + std::string(argv[2]) + std::string(argv[3]) + ".cut");
  DeltaEFile << Lower << " " << Upper;
  DeltaEFile.close();
  std::cout << "deltaE cuts saved\n";
  return 0;
}
