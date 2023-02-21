// Martin Duy Tat 15th February 2023

#include<vector>
#include<numeric>
#include<string>
#include<iostream>
#include<iomanip>
#include"uncertainties/ureal.hpp"
#include"uncertainties/io.hpp"
#include"Settings.h"
#include"Utilities.h"
#include"FlavourTags/KiCombiner.h"

KiCombiner::KiCombiner(const Settings &settings):
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")),
  m_FlavourTags(LoadFlavourTags(settings)) {
}

std::vector<uncertainties::udouble> KiCombiner::GetKiWithUncertainties(
  const std::vector<double> &ci,
  const std::vector<double> &si) const {
  std::vector<uncertainties::udouble> Ki(2*m_NumberBins);
  for(const auto &Tag : m_FlavourTags) {
    std::size_t Index = 0;
    for(int Bin = -m_NumberBins; Bin <= static_cast<int>(m_NumberBins); Bin++) {
      if(Bin == 0) {
	continue;
      }
      Ki[Index] += Tag.GetKi(Bin, ci, si);
      Index++;
    }
  }
  const uncertainties::udouble Sum = std::accumulate(Ki.begin(),
						     Ki.end(),
						     uncertainties::udouble(0.0, 0.0));
  std::transform(Ki.begin(), Ki.end(), Ki.begin(),
		 [=] (const auto &a) { return a/Sum; });
  return Ki;
}

std::vector<double> KiCombiner::GetKi(const std::vector<double> &ci,
				      const std::vector<double> &si) const {
  const auto Ki_unc = GetKiWithUncertainties(ci, si);
  std::vector<double> Ki;
  for(const auto Ki_u : Ki_unc) {
    Ki.push_back(uncertainties::nom(Ki_u));
  }
  return Ki;
}

void KiCombiner::PrintKi(const std::vector<double> &ci,
			 const std::vector<double> &si) const {
  std::cout << std::left << std::setw(10) << "Bin";
  std::cout << std::left << std::setw(20) << "Ki";
  std::cout << std::left << std::setw(20) << "Kbari";
  std::cout << "\n";
  const auto Ki = GetKiWithUncertainties(ci, si);
  const std::size_t Size = ci.size();
  for(std::size_t Bin = 1; Bin <= Size; Bin++) {
    std::cout << std::left << std::setw(10) << Bin;
    std::cout << std::left << std::setw(20) << Ki[Bin - 1 + Size];
    std::cout << std::left << std::setw(20) << Ki[-static_cast<int>(Bin) + Size];
    std::cout << "\n";
  }
  for(const auto &Tag : m_FlavourTags) {
    Tag.PrintDCS(ci, si);
  }
}

std::vector<FlavourTagYield>
KiCombiner::LoadFlavourTags(const Settings &settings) const {
  std::vector<std::string> TagModes =
    Utilities::ConvertStringToVector(settings.get("FlavourTagModes"));
  std::vector<FlavourTagYield> FlavourTags;
  for(const auto Tag : TagModes) {
    FlavourTags.push_back(FlavourTagYield(Tag, settings));
  }
  return FlavourTags;
}
