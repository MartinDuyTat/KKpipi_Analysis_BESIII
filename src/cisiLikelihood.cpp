// Martin Duy Tat 13th October 2022

#include<vector>
#include<algorithm>
#include<string>
#include<iterator>
#include<numeric>
#include<utility>
#include<iostream>
#include<iomanip>
#include"Settings.h"
#include"cisiLikelihood.h"
#include"BinnedDTData.h"
#include"Utilities.h"
#include"cisiFitterParameters.h"

cisiLikelihood::cisiLikelihood(const Settings &settings):
  m_TagData(SetupTags(settings)) {
}

double cisiLikelihood::CalculateLogLikelihood(
  const cisiFitterParameters &Parameters) const {
  auto LikelihoodAdder = [&] (double a, const BinnedDTData &b) {
    return a + b.GetLogLikelihood(Parameters);
  };
  return std::accumulate(m_TagData.begin(), m_TagData.end(), 0.0, LikelihoodAdder);
}

void cisiLikelihood::LoadToyDataset(int ToyNumber) const {
  for(const auto &TagData : m_TagData) {
    TagData.LoadToyDataset(ToyNumber);
  }
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

void cisiLikelihood::PrintComparison(const cisiFitterParameters &Parameters) const {
  std::for_each(m_TagData.begin(),
		m_TagData.end(),
		[&] (const auto &a) {
		  a.PrintComparison(Parameters); });
  PrintFinalKi(Parameters.m_Ri);
}			     

void cisiLikelihood::PrintFinalKi(const std::vector<double> &Ri) const {
  std::vector<double> Ki, Kbari;
  Utilities::ConvertRiToKi(Ri, Ki, Kbari);
  std::cout << std::left << std::setw(10) << "Ki";
  std::cout << std::left << std::setw(10) << "Kbari" << "\n";
  for(std::size_t i = 0; i < Ki.size(); i++) {
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << Ki[i];
    std::cout << std::left << std::setw(10);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << Kbari[i] << "\n";
  }
}

void cisiLikelihood::SavePredictedBinYields(
  std::ofstream &File,
  const cisiFitterParameters &Parameters) const {
  std::for_each(m_TagData.begin(),
		m_TagData.end(),
		[&] (const auto &a) {
		  a.SavePredictedBinYields(File, Parameters); });
}
