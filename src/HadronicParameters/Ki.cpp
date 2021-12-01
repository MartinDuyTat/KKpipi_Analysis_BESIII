// Martin Duy Tat 1st December 2021

#include<stdexcept>
#include"Settings.h"
#include"HadronicParameters/Ki.h"

Ki::Ki(const Settings &settings): m_Ki(settings.getI("NumberBins")), m_Kbari(settings.getI("NumberBins")) {
  for(int i = 1; i <= settings.getI("NumberBins"); i++) {
    m_Ki[i - 1] = settings["Ki"].getD("Model_K_p" + std::to_string(i));
    m_Kbari[i - 1] = settings["Ki"].getD("Model_K_m" + std::to_string(i));
  }
}

double Ki::Get_Ki(int Bin) const {
  if(Bin == 0) {
    throw std::out_of_range("Bin number cannot be 0");
  }
  return Bin > 0 ? m_Ki[Bin - 1] : m_Kbari[-Bin - 1];
}
