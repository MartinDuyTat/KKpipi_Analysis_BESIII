// Martin Duy Tat 3rd October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixTSym.h"
#include"RawBinnedCPTagYields.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"
#include"Utilities.h"

RawBinnedCPTagYields::RawBinnedCPTagYields(const std::string &Tag,
					   const Settings &settings):
  RawBinnedDTYields(ParseYields(Tag, settings),
		    LoadCorrelationMatrix(Tag, settings)) {
}

std::vector<AsymmetricUncertainty> 
RawBinnedCPTagYields::ParseYields(const std::string &Tag,
				  const Settings &settings) const {
  std::vector<AsymmetricUncertainty> Yields;
  const std::string DTYieldFilename =
    Utilities::ReplaceString(settings.get("DT_Yield"), "TAG", Tag);
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

TMatrixTSym<double>
RawBinnedCPTagYields::LoadCorrelationMatrix(const std::string &Tag,
					    const Settings &settings) const {
  const std::string Filename =
    Utilities::ReplaceString(settings.get("RawYields_CorrMatrix"), "TAG", Tag);
  TFile File(Filename.c_str());
  TMatrixTSym<double> *CorrelationMatrix = nullptr;
  File.GetObject("CorrelationMatrix", CorrelationMatrix);
  File.Close();
  return *CorrelationMatrix;
}
