// Martin Duy Tat 25th November 2021

#include<string>
#include<vector>
#include<stdexcept>
#include<numeric>
#include<algorithm>
#include<utility>
#include"TMath.h"
#include"Category.h"
#include"Settings.h"

Category::Category(const Settings &settings): m_TagMode(settings.get("Mode")),
					      m_SignalBins(settings["BinningScheme"].getI("NumberBins")),
					      m_Inclusive(settings.contains("Inclusive_fit") && settings.getB("Inclusive_fit")),
					      m_CategoryVar(("DoubleTag_" + m_TagMode + "_Categories").c_str(), "") {
  if(settings.get("TagType") != "DT") {
    throw std::invalid_argument("Cannot do double tag fit with single tags!");
  }
  if(m_Inclusive) {
    m_SignalBins = 0;
  }
  if(m_TagMode == "Kpi" || m_TagMode == "Kpipi0" || m_TagMode == "Kpipipi" || m_TagMode == "KeNu") {
    m_Type = "Flavour";
    m_TagBins = 0;
  } else if(m_TagMode == "KSpipi" || m_TagMode == "KLpipi" || m_TagMode == "KSpipiPartReco") {
    m_Type = "SCMB";
    m_TagBins = 8;
  } else if(m_TagMode == "KSKK") {
    m_Type = "SCMB";
    m_TagBins = 2;
  } else if(m_TagMode == "KKpipi") {
    m_Type = "SCMB";
    m_TagBins = 8;
  } else {
    m_Type = "CP";
    m_TagBins = 0;
  }
  for(const auto &Category : GetCategories()) {
    m_CategoryVar.defineType(Category.c_str());
  }
}

void Category::CheckValidBins(int SignalBin, int TagBin) const {
  // For inclusive fits, both bins must be zero
  if(m_Inclusive) {
    if(SignalBin != 0) {
      throw std::out_of_range("Inclusive fit cannot have non-zero siganl bin!\n");
    } else {
      return;
    }
  }
  // For binned fits it's much more complicated...
  // Signal KKpipi bin must be +-1, +-2, ..., +- 8
  if(SignalBin == 0 || TMath::Abs(SignalBin) > m_SignalBins) {
    throw std::out_of_range("Signal bin number " + std::to_string(SignalBin) + " does not exist!");
  }
  // For flavour or CP tags, there is no tag side binning
  if((m_Type == "Flavour" || m_Type == "CP") && TagBin != 0) {
    throw std::out_of_range("Tag bin number " + std::to_string(TagBin) + " does not exist for " + m_TagMode + " tag mode!");
  }
  // For SCMB tags, tag side bins must be positive and within the correct range
  if((m_Type == "KSpipi" || m_Type == "KSKK" || m_Type == "KKpipi") && (TagBin <= 0 || TagBin > m_TagBins)) {
    throw std::out_of_range("Tag bin number " + std::to_string(TagBin) + " does not exist for " + m_TagMode + " tag mode!");
  }
}

std::string Category::GetCategory(int SignalBin, int TagBin) const {
  // Make sure bins are valid
  CheckValidBins(SignalBin, TagBin);
  // Start building up the category string
  std::string CategoryString("DoubleTag_");
  CategoryString += m_Type;
  CategoryString += "_KKpipi_vs_" + m_TagMode;
  // Signal bin number
  CategoryString += "_SignalBin";
  if(!m_Inclusive && (m_Type == "Flavour" || m_Type == "SCMB")) {
    CategoryString += SignalBin > 0 ? "P" : "M";
  }
  CategoryString += std::to_string(TMath::Abs(SignalBin));
  // Tag bin number
  if(m_Type != "CP") {
    CategoryString += "_TagBin" + std::to_string(TagBin);
  }
  return CategoryString;
}

std::string Category::operator ()(int SignalBin, int TagBin) const {
  return GetCategory(SignalBin, TagBin);
}

std::vector<std::pair<int, int>> Category::GetBinCombinations() const {
  std::vector<int> SignalBins(m_SignalBins), TagBins(m_TagBins);
  // Signal side binning
  if(m_Inclusive) {
    // For inclusive fits we only have a single category
    SignalBins.push_back(0);
  } else {
    // For binned fits it's much more complicated...
    // We need at least 8 bins on the signal KKpipi side
    std::iota(SignalBins.begin(), SignalBins.end(), 1);
    // If tag mode is not KKpipi, we also need the conjugate bins
    if(m_TagMode != "KKpipi") {
      // Double the number of bins
      SignalBins.resize(2*m_SignalBins);
      // Fill with negative bins
      std::iota(SignalBins.begin() + m_SignalBins, SignalBins.end(), -m_SignalBins);
      // Swap the ordering so that we have bins from -8, -7, ..., -1, 1, 2, ..., 8
      std::rotate(SignalBins.begin(), SignalBins.begin() + m_SignalBins, SignalBins.end());
    }
  }
  // Tag side binning
  if(m_Type == "Flavour" || m_Type == "CP") {
    // For CP and flavour tags, there is no binning on the tag side
    TagBins.push_back(0);
  } else {
    // For SCMB tags we need binning on the tag side
    std::iota(TagBins.begin(), TagBins.end(), 1);
  }
  // Combine all bin combinations in the correct order
  std::vector<std::pair<int, int>> BinCombinations;
  // Loop over tag side bins
  for(int j : TagBins) {
    // Loop over signal side bins
    for(int i : SignalBins) {
      // For KKpipi vs KKpipi tags, bin ij is equivalent to ji
      if(m_TagMode == "KKpipi" && i < j) {
	continue;
      }
      BinCombinations.push_back({i, j});
    }
  }
  return BinCombinations;
}

std::vector<std::string> Category::GetCategories() const {
  std::vector<std::string> CategoryStrings;
  for(const auto &BinCombination : GetBinCombinations()) {
    CategoryStrings.push_back(GetCategory(BinCombination.first, BinCombination.second));
  }
  return CategoryStrings;
}

RooCategory* Category::GetCategoryVariable() {
  return &m_CategoryVar;
}

int Category::GetBinNumber(const std::string &category, const std::string &SignalTag) const {
  if(SignalTag != "Signal" && SignalTag != "Tag") {
    return 0;
  }
  // Find position of the bin number in the category string
  auto pos = category.find(SignalTag + "Bin") + SignalTag.length() + 3;
  // Check if next character is P for plus, M for minus or just a bin number
  int Sign;
  if(category[pos] == 'P') {
    Sign = +1;
  } else if(category[pos] == 'M') {
    Sign = -1;
  } else {
    // If there is no P or M, then this tag doesn't have conjugate bins
    return static_cast<int>(category[pos] - '0');
  }
  // Convert char number to int with some char trickery
  int Number = category[pos + 1] - '0';
  return Sign*Number;
}

int Category::GetSignalBinNumber(const std::string &category) const {
  // Inclusive fits only have a single bin
  if(m_Inclusive) {
    return 0;
  } else {
    return GetBinNumber(category, "Signal");
  }
}

int Category::GetTagBinNumber(const std::string &category) const {
  return GetBinNumber(category, "Tag");
}
