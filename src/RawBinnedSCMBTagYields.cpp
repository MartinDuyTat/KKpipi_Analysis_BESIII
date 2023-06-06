// Martin Duy Tat 3rd October 2022

#include<vector>
#include<string>
#include"RawBinnedSCMBTagYields.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedSCMBTagYields::RawBinnedSCMBTagYields(const std::string &Tag,
					       const Settings &settings,
					       int ToyNumber):
  RawBinnedDTYields(ParseYields(Tag, settings, ToyNumber),
		    LoadCorrelationMatrix(Tag, settings)) {
}

std::vector<AsymmetricUncertainty> 
RawBinnedSCMBTagYields::ParseYields(const std::string &Tag,
				    const Settings &settings,
				    int ToyNumber) const {
  std::vector<AsymmetricUncertainty> Yields;
  const auto DTYieldFilename = GetFilename(Tag, settings, ToyNumber);
  const auto ParsedDTYields = Utilities::ParseFile(DTYieldFilename);
  const std::string YieldNamePrefix("DoubleTag_SCMB_KKpipi_vs_" + Tag);
  std::size_t NumberBins = settings["BinningScheme"].getI("NumberBins");
  for(std::size_t TagBin = 1; TagBin <= 8; TagBin++) {
    for(int SignalBin = -NumberBins;
	SignalBin <= static_cast<int>(NumberBins);
	SignalBin++) {
      if(SignalBin == 0) {
	continue;
      }
      const std::string YieldName(YieldNamePrefix +
				  "_SignalBin" +
				  (SignalBin > 0 ? "P" : "M") +
				  std::to_string(TMath::Abs(SignalBin)) +
				  "_TagBin" + 
				  std::to_string(TagBin) +
				  "_SignalYield");
      const double Yield = ParsedDTYields.at(YieldName);
      const double PlusError = ParsedDTYields.at(YieldName + "_high_err");
      const double MinusError = -ParsedDTYields.at(YieldName + "_low_err");
      const double SymmetricError = ParsedDTYields.at(YieldName + "_err");
      Yields.push_back({Yield, PlusError, MinusError, SymmetricError});
    }
  }
  return Yields;
}
