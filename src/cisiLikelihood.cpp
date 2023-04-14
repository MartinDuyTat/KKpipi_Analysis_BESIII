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
  m_TagData(SetupTags(settings)) {
}

double cisiLikelihood::CalculateLogLikelihood(double BF_KKpipi,
					      const std::vector<double> &ci,
					      const std::vector<double> &si,
					      const std::vector<double> &Ki,
					      const std::vector<double> &Kbari,
					      double DeltaKpi) const {
  auto LikelihoodAdder = [&] (double a, const BinnedDTData &b) {
    return a + b.GetLogLikelihood(BF_KKpipi, ci, si, Ki, Kbari, DeltaKpi);
  };
  return std::accumulate(m_TagData.begin(), m_TagData.end(), 0.0, LikelihoodAdder);
}

void cisiLikelihood::GenerateToy(double BF_KKpipi,
				 const std::vector<double> &ci,
				 const std::vector<double> &si,
				 const std::vector<double> &Ki,
				 const std::vector<double> &Kbari,
				 double DeltaKpi,
				 std::size_t StatsMultiplier) const {
  for(const auto &TagData : m_TagData) {
    TagData.GenerateToyYields(BF_KKpipi,
			      ci, si,
			      Ki, Kbari,
			      DeltaKpi,
			      StatsMultiplier);
  }
}

double cisiLikelihood::CalculateToyLogLikelihood(
  double BF_KKpipi,
  const std::vector<double> &ci,
  const std::vector<double> &si,
  const std::vector<double> &Ki,
  const std::vector<double> &Kbari,
  double DeltaKpi) const {
  auto LikelihoodAdder = [&] (double a, const BinnedDTData &b) {
    return a + b.GetToyLogLikelihood(BF_KKpipi, ci, si, Ki, Kbari, DeltaKpi);
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
		   return BinnedDTData(Tag, settings);
		 });
  return TagData;
}

void cisiLikelihood::PrintComparison(double BF_KKpipi,
				     const std::vector<double> &ci,
				     const std::vector<double> &si,
				     const std::vector<double> &Ki,
				     const std::vector<double> &Kbari,
				     double DeltaKpi) const {
  std::for_each(m_TagData.begin(),
		m_TagData.end(),
		[&] (const auto &a) {
		  a.PrintComparison(BF_KKpipi, ci, si, Ki, Kbari, DeltaKpi); });
}			     

/*std::pair<std::vector<double>, std::vector<double>>
cisiLikelihood::GetKi(const Settings &settings) const {
  int NumberBins = settings["BinningScheme"].getI("NumberBins");
  std::vector<double> Ki, Kbari;
  for(int Bin = 1; Bin <= NumberBins; Bin++) {
    const std::string KiName = "KKpipi_Ki_SignalBinP" + std::to_string(Bin);
    const std::string KbariName = "KKpipi_Ki_SignalBinM" + std::to_string(Bin);
    Ki.push_back(settings["FractionalYields"].getD(KiName));
    Kbari.push_back(settings["FractionalYields"].getD(KbariName));
  }
  const double Sum = std::accumulate(Ki.begin(), Ki.end(), 0.0)
                   + std::accumulate(Kbari.begin(), Kbari.end(), 0.0);
  std::transform(Ki.begin(), Ki.end(), Ki.begin(),
		 [=] (double a) {return a/Sum;});
  std::transform(Kbari.begin(), Kbari.end(), Kbari.begin(),
		 [=] (double a) {return a/Sum;});
  return std::make_pair(Ki, Kbari);
}*/
