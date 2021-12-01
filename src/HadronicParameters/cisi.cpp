// Martin Duy Tat 1st December 2021

#include<stdexcept>
#include"Settings.h"
#include"HadronicParameters/cisi.h"

cisi::cisi(const Settings &settings): m_ci(settings.getI("NumberBins")), m_si(settings.getI("NumberBins")) {
  for(int i = 1; i <= settings.getI("NumberBins"); i++) {
    m_ci[i - 1] = settings["cisi"].getD("Model_c" + std::to_string(i));
    m_si[i - 1] = settings["cisi"].getD("Model_s" + std::to_string(i));
  }
}

double cisi::Get_ci(int Bin) const {
  if(Bin == 0) {
    throw std::out_of_range("Bin number cannot be 0");
  }
  return Bin > 0 ? m_ci[Bin - 1] : m_ci[-Bin - 1];
}

double cisi::Get_si(int Bin) const {
  if(Bin == 0) {
    throw std::out_of_range("Bin number cannot be 0");
  }
  return Bin > 0 ? m_si[Bin - 1] : -m_si[-Bin - 1];
}
