// Martin Duy Tat 3rd October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixTSym.h"
#include"RawBinnedFlavourTagYields.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedFlavourTagYields::RawBinnedFlavourTagYields(const std::string &Tag,
					   const Settings &settings):
  RawBinnedDTYields(ParseYields(Tag, settings),
		    LoadCorrelationMatrix(Tag, settings)) {
}

std::vector<AsymmetricUncertainty> 
RawBinnedFlavourTagYields::ParseYields(const std::string &Tag,
				  const Settings &settings) const {
  std::vector<AsymmetricUncertainty> Yields;
  const std::string DTYieldFilename =
    Utilities::ReplaceString(settings.get("DT_Yield"), "TAG", Tag);
  const auto ParsedDTYields = Utilities::ParseFile(DTYieldFilename);
  const std::string YieldNamePrefix("DoubleTag_Flavour_KKpipi_vs_" + Tag);
  std::size_t NumberBins = settings["BinningScheme"].getI("NumberBins");
  for(int Bin = -NumberBins; Bin <= static_cast<int>(NumberBins); Bin++) {
    if(Bin == 0) {
      continue;
    }
    const std::string YieldName(YieldNamePrefix +
				"_SignalBin" +
				(Bin > 0 ? "P" : "M") + 
				std::to_string(TMath::Abs(Bin)) +
				"_TagBin0" +
				"_SignalYield");
    const double Yield = ParsedDTYields.at(YieldName);
    const double PlusError = ParsedDTYields.at(YieldName + "_high_err");
    const double MinusError = -ParsedDTYields.at(YieldName + "_low_err");
    const double SymmetricError = ParsedDTYields.at(YieldName + "_err");
    Yields.push_back({Yield, PlusError, MinusError, SymmetricError});
  }
  return Yields;
}

TMatrixTSym<double>
RawBinnedFlavourTagYields::LoadCorrelationMatrix(const std::string &Tag,
					    const Settings &settings) const {
  const std::string Filename =
    Utilities::ReplaceString(settings.get("RawYields_CorrMatrix"), "TAG", Tag);
  TFile File(Filename.c_str());
  TMatrixTSym<double> *CorrelationMatrix = nullptr;
  File.GetObject("CorrelationMatrix", CorrelationMatrix);
  File.Close();
  return *CorrelationMatrix;
}
