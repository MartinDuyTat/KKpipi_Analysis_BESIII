// Martin Duy Tat 13th October 2022

#include<vector>
#include<algorithm>
#include<string>
#include<iterator>
#include<numeric>
#include<utility>
#include"Settings.h"
#include"cisiLikelihood.h"
#include"BinnedDTData.h"
#include"Utilities.h"

cisiLikelihood::cisiLikelihood(const Settings &settings):
  m_TagData(SetupTags(settings)),
  m_Ki(GetKi(settings).first),
  m_Kbari(GetKi(settings).second) {
}

double cisiLikelihood::CalculateLogLikelihood(const std::vector<double> &ci,
					      const std::vector<double> &si) const {
  auto LikelihoodAdder = [&] (double a, const BinnedDTData &b) {
    return a + b.GetLogLikelihood(ci, si);
  };
  return std::accumulate(m_TagData.begin(), m_TagData.end(), 0.0, LikelihoodAdder);
}

std::vector<BinnedDTData> cisiLikelihood::SetupTags(const Settings &settings) const {
  const std::vector<std::string> TagModes = 
    Utilities::ConvertStringToVector(settings.get("TagModes"));
  std::vector<BinnedDTData> TagData;
  std::transform(TagModes.begin(),
		 TagModes.end(),
		 std::back_inserter(TagData),
		 [&] (const auto &Tag) {
		   return BinnedDTData(Tag, m_Ki, m_Kbari, settings);
		 });
  return TagData;
}

void cisiLikelihood::PrintComparison(const std::vector<double> &ci,
				     const std::vector<double> &si) const {
  std::for_each(m_TagData.begin(),
		m_TagData.end(),
		[&] (const auto &a) { a.PrintComparison(ci, si); });
}
			     

std::pair<std::vector<double>, std::vector<double>>
cisiLikelihood::GetKi(const Settings &settings) const {
  int NumberBins = settings["BinningScheme"].getI("NumberBins");
  std::vector<double> Ki, Kbari;
  for(int Bin = 1; Bin <= NumberBins; Bin++) {
    const std::string KiName = "KKpipi_Ki_SignalBinP" + std::to_string(Bin);
    const std::string KbariName = "KKpipi_Ki_SignalBinM" + std::to_string(Bin);
    Ki.push_back(settings["FractionalYields"].getD(KiName));
    Kbari.push_back(settings["FractionalYields"].getD(KbariName));
  }
  return std::make_pair(Ki, Kbari);
}
