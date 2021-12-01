// Martin Duy Tat 1st December 2021

#include<algorithm>
#include<sstream>
#include"HadronicParameters/DCS_Parameters.h"
#include"Settings.h"
#include"uncertainties/impl.hpp"
#include"uncertainties/ureal.hpp"
#include"uncertainties/ureals.hpp"

using uncertainties::udouble;
using uncertainties::ureals;

DCS_Parameters::DCS_Parameters(const Settings &settings): m_DCS(3), m_DCS_Cov(9) {
  // List of parameters
  std::vector<std::string> Parameters{"rD", "R", "deltaD"};
  // Parse central values
  for(int i = 0; i < 3; i++) {
    m_DCS[i] = settings.getD(Parameters[i]);
  }
  // Parse correlation matrix
  std::string CovString = settings.get("CorrelationMatrix");
  std::replace(CovString.begin(), CovString.end(), ',', ' ');
  std::stringstream ss(CovString);
  for(int i = 0; i < 9; i++) {
    double Value;
    ss >> Value;
    m_DCS_Cov[i] = Value;
  }
  // Convert correlation matrix to covariance matrix
  for(int i = 0 ; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m_DCS_Cov[3*i + j] *= settings.getD(Parameters[i] + "_err")*settings.getD(Parameters[j] + "_err");
    }
  }
}

std::vector<udouble> DCS_Parameters::GetDCSParameters() const {
  return ureals<std::vector<udouble>>(m_DCS, m_DCS_Cov);
}
