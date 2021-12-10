// Martin Duy Tat 4th April 2021

#include<string>
#include"TCut.h"
#include"TruthMatchingCuts.h"
#include"CutsFromFile.h"

TruthMatchingCuts::TruthMatchingCuts(const std::string &SignalMode,
				     const std::string &TagMode): m_SignalMode(SignalMode), m_TagMode(TagMode) {
}

TCut TruthMatchingCuts::CheckEmptyCut(TCut Cut, const std::string &TagSide) const {
  if(Cut == TCut()) {
    Cut = TCut((TagSide + "IsSameDMother == 1 && " + TagSide + "PIDTrue == 1").c_str());
  }
  return Cut;
}

TCut TruthMatchingCuts::GetTruthMatchingCuts() const {
  if(m_SignalMode == "") {
    CutsFromFile cutsFromFile(std::string(TRUTH_MATCHING_CUTS_DIR) + m_TagMode + ".cut", "");
    return CheckEmptyCut(cutsFromFile.GetCuts());
  } else {
    CutsFromFile cutsFromFileSignal(std::string(TRUTH_MATCHING_CUTS_DIR) + m_SignalMode + ".cut", "Signal");
    CutsFromFile cutsFromFileTag(std::string(TRUTH_MATCHING_CUTS_DIR) + m_TagMode + ".cut", "Tag");
    return CheckEmptyCut(cutsFromFileSignal.GetCuts(), "Signal") && CheckEmptyCut(cutsFromFileTag.GetCuts(), "Tag");
  }
}
