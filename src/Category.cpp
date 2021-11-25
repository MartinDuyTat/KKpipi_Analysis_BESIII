// Martin Duy Tat 25th November 2021

#include<string>
#include<vector>
#include<stdexcept>
#include<numeric>
#include"TMath.h"
#include"Category.h"
#include"Settings.h"

#include<iostream> //remove

Category::Category(const Settings &settings): m_TagMode(settings.get("Mode")),
					      m_SignalBins(4),
					      m_CategoryVar(("DoubleTag_" + m_TagMode + "_Categories").c_str(), "") {
  if(settings.get("TagType") != "DT") {
    throw std::invalid_argument("Cannot do binned fit with single tags!");
  }
  if(m_TagMode == "Kpi" || m_TagMode == "Kpipi0" || m_TagMode == "Kpipipi" || m_TagMode == "KeNu") {
    m_Type = "Flavour";
    m_TagBins = 0;
  } else if(m_TagMode == "KSpipi") {
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

std::string Category::GetCategory(int SignalBin, int TagBin) const {
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
  std::string CategoryString("DoubleTag_");
  CategoryString += m_Type + "_";
  CategoryString += "KKpipi_vs_" + m_TagMode + "_";
  CategoryString += "SignalBin";
  if(m_Type == "Flavour" || m_Type == "SCMB") {
    CategoryString += SignalBin > 0 ? "P" : "M";
  }
  CategoryString += std::to_string(TMath::Abs(SignalBin));
  if(m_Type != "CP") {
    CategoryString += "_TagBin" + std::to_string(TagBin);
  }
  return CategoryString;
}

std::string Category::operator ()(int SignalBin, int TagBin) const {
  return GetCategory(SignalBin, TagBin);
}

std::vector<std::string> Category::GetCategories() const {
  std::vector<std::string> CategoryStrings;
  std::vector<int> SignalBins(8), TagBins(m_TagBins);
  // We need at least 8 bins on the signal KKpipi side
  std::iota(SignalBins.begin(), SignalBins.end(), 1);
  // If tag mode is not KKpipi, we also need the conjugate bins
  if(m_TagMode != "KKpipi") {
    for(int i = 1; i <= 8; i++) {
      SignalBins.push_back(-i);
    }
  }
  if(m_Type == "Flavour" || m_Type == "CP") {
    // For CP and flavour tags, there is no binning on the tag side
    TagBins.push_back(0);
  } else {
    // For SCMB tags we need binning on the tag side
    std::iota(TagBins.begin(), TagBins.end(), 1);
  }
  // Loop over signal side bins
  for(int i : SignalBins) {
    // Loop over tag side bins
    for(int j : TagBins) {
      // For KKpipi vs KKpipi tags, bin ij is equivalent to ji
      if(m_TagMode == "KKpipi" && i < j) {
	continue;
      }
      CategoryStrings.push_back(GetCategory(i, j));
      // For cross check, remove
      std::cout << CategoryStrings.back() << "\n";
    }
  }
  return CategoryStrings;
}

RooCategory Category::GetCategoryVariable() const {
  return m_CategoryVar;
}
