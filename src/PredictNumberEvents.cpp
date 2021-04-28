// Martin Duy Tat 28th April 2021

#include<vector>
#include"TMath.h"
#include"PredictNumberEvents.h"
#include"DDecayParameters.h"

PredictNumberEvents::PredictNumberEvents(const DDecayParameters &DDParameters, double DTFlavourYield, double STFlavourYield, double STTagYield, int CP): m_K(DDParameters.GetK()), m_Kbar(DDParameters.GetKbar()), m_ci(DDParameters.Getc()), m_si(DDParameters.Gets()), m_NBins(static_cast<int>(m_K.size())), m_DTFlavourYield(DTFlavourYield), m_STFlavourYield(STFlavourYield), m_STTagYield(STTagYield), m_CP(CP) {
}

double PredictNumberEvents::GetCPEvenTagYield(int i) const {
  if(TMath::Abs(i) > m_NBins || i == 0) {
    return 0;
  }
  int BinIndex = TMath::Abs(i) - 1;
  return m_DTFlavourYield*m_STTagYield*(m_K[BinIndex] + m_Kbar[BinIndex] - 2*m_ci[BinIndex]*TMath::Sqrt(m_K[BinIndex]*m_Kbar[BinIndex]))/(2.0*m_STFlavourYield);
}

double PredictNumberEvents::GetCPOddTagYield(int i) const {
  if(TMath::Abs(i) > m_NBins || i == 0) {
    return 0;
  }
  int BinIndex = TMath::Abs(i) - 1;
  return m_DTFlavourYield*m_STTagYield*(m_K[BinIndex] + m_Kbar[BinIndex] + 2*m_ci[BinIndex]*TMath::Sqrt(m_K[BinIndex]*m_Kbar[BinIndex]))/(2.0*m_STFlavourYield);
}

double PredictNumberEvents::GetFlavourTagYield(int i) const {
  if(TMath::Abs(i) > m_NBins || i == 0) {
    return 0;
  }
  int BinIndex = TMath::Abs(i) - 1;
  if(i > 0) {
    return m_DTFlavourYield*m_K[BinIndex];
  } else {
    return m_DTFlavourYield*m_Kbar[BinIndex];
  }
}

double PredictNumberEvents::operator()(int i) const {
  switch(m_CP) {
  case +1:
    return GetCPEvenTagYield(i);
    break;
  case -1:
    return GetCPOddTagYield(i);
    break;
  case 0:
    return GetFlavourTagYield(i);
    break;
  default:
    return 0;
  }
}

int PredictNumberEvents::GetNBins() const {
  return m_NBins;
}
