// Martin Duy Tat 4th April 2021

#include<string>
#include"TCut.h"
#include"TruthMatchingCuts.h"

TruthMatchingCuts::TruthMatchingCuts(const std::string &TagType): m_TagType(TagType) {
}

TCut TruthMatchingCuts::GetTruthMatchingCuts() const {
  if(m_TagType == "ST") {
    return TCut("IsSameDMother == 1 && PIDTrue == 1");
  } else if(m_TagType == "DT") {
    return TCut("SignalIsSameDMother == 1 && SignalPIDTrue == 1 && TagIsSameDMother == 1 && PIDTrue == 1");
  } else {
    return TCut();
  }
}
