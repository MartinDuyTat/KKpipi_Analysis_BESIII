// Martin Duy Tat 3rd October 2022

#include<vector>
#include<string>
#include"RawBinnedCPTagYields.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedCPTagYields::RawBinnedCPTagYields(const std::string &Tag,
					   const Settings &settings,
					   int ToyNumber):
  RawBinnedDTYields(ParseYields(Tag, settings, ToyNumber),
		    LoadCorrelationMatrix(Tag, settings)) {
}

std::vector<AsymmetricUncertainty> 
RawBinnedCPTagYields::ParseYields(const std::string &Tag,
				  const Settings &settings,
				  int ToyNumber) const {
  std::vector<AsymmetricUncertainty> Yields;
  std::string DTYieldFilename;
  if(ToyNumber > 0) {
    DTYieldFilename =
      Utilities::ReplaceString(settings.get("DT_ToyYieldDir"), "TAG", Tag);
    DTYieldFilename += "/Toy" + std::to_string(ToyNumber) + ".txt";
  } else {
    DTYieldFilename =
      Utilities::ReplaceString(settings.get("DT_Yield"), "TAG", Tag);
  }
  const auto ParsedDTYields = Utilities::ParseFile(DTYieldFilename);
  const std::string YieldNamePrefix("DoubleTag_CP_KKpipi_vs_" + Tag);
  std::size_t NumberBins = settings["BinningScheme"].getI("NumberBins");
  for(std::size_t Bin = 1; Bin <= NumberBins; Bin++) {
    const std::string YieldName(YieldNamePrefix +
				"_SignalBin" +
				std::to_string(Bin) +
				"_SignalYield");
    const double Yield = ParsedDTYields.at(YieldName);
    const double PlusError = ParsedDTYields.at(YieldName + "_high_err");
    const double MinusError = -ParsedDTYields.at(YieldName + "_low_err");
    const double SymmetricError = ParsedDTYields.at(YieldName + "_err");
    Yields.push_back({Yield, PlusError, MinusError, SymmetricError});
  }
  return Yields;
}
