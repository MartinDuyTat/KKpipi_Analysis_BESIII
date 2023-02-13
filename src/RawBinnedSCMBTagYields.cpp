// Martin Duy Tat 3rd October 2022

#include<vector>
#include<string>
#include"TFile.h"
#include"TMatrixTSym.h"
#include"RawBinnedSCMBTagYields.h"
#include"RawBinnedDTYields.h"
#include"Settings.h"

RawBinnedSCMBTagYields::RawBinnedSCMBTagYields(const std::string &Tag,
					   const Settings &settings):
  RawBinnedDTYields(ParseYields(Tag, settings),
		    LoadCorrelationMatrix(Tag, settings)) {
}

std::vector<AsymmetricUncertainty> 
RawBinnedSCMBTagYields::ParseYields(const std::string &Tag,
				  const Settings &settings) const {
  std::vector<AsymmetricUncertainty> Yields;
  const std::string SettingsName(Tag + "_DT_Yield");
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
      const double Yield = settings[SettingsName].getD(YieldName);
      const double PlusError = settings[SettingsName].getD(YieldName + "_high_err");
      const double MinusError = -settings[SettingsName].getD(YieldName + "_low_err");
      Yields.push_back({Yield, PlusError, MinusError});
    }
  }
  return Yields;
}

TMatrixTSym<double>
RawBinnedSCMBTagYields::LoadCorrelationMatrix(const std::string &Tag,
					    const Settings &settings) const {
  const std::string Filename = settings.get(Tag + "_RawYields_CorrelationMatrix");
  TFile File(Filename.c_str());
  TMatrixTSym<double> *CorrelationMatrix = nullptr;
  File.GetObject("CorrelationMatrix", CorrelationMatrix);
  File.Close();
  return *CorrelationMatrix;
}
