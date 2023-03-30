// Martin Duy Tat 15th February 2023

#include<stdexcept>
#include<string>
#include<vector>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<iomanip>
#include"uncertainties/ureal.hpp"
#include"uncertainties/ureals.hpp"
#include"TMatrixT.h"
#include"TMatrixTSym.h"
#include"TMath.h"
#include"TFile.h"
#include"Settings.h"
#include"FlavourTags/FlavourTagYield.h"

FlavourTagYield::FlavourTagYield(const std::string &TagMode,
				 const Settings &settings):
  m_TagMode(TagMode),
  m_IncludeDCSCorrections(settings.getB("DCS_Corrections") && TagMode != "KeNu"),
  m_DTFlavourYields(LoadFlavourYields(TagMode, settings)),
  m_rD(LoadCharmParameter(TagMode, settings, "rD")),
  m_R(LoadCharmParameter(TagMode, settings, "R")),
  m_DeltaD(LoadCharmParameter(TagMode, settings, "deltaD")),
  m_SinDeltaD(TMath::Sin(TMath::Pi()*m_DeltaD.first/180.0)),
  m_CosDeltaD(TMath::Cos(TMath::Pi()*m_DeltaD.first/180.0)),
  m_CharmCovMatrix(LoadCharmCovMatrix(TagMode, settings)),
  m_RawSingleTagYield(GetSTYield(TagMode, settings)),
  m_SingleTagEfficiency(GetSTEfficiency(TagMode, settings)) {
}   

uncertainties::udouble FlavourTagYield::GetKi(
  int Bin,
  const std::vector<double> &ci,
  const std::vector<double> &si) const {
  const double DCSCorrection = 
    m_IncludeDCSCorrections ? CalculateDCSCorrection(Bin, ci, si) : 1.0;
  const uncertainties::udouble Yield = GetDTFlavourYield(Bin);
  const double SingleTagYield =
    m_RawSingleTagYield.first/m_SingleTagEfficiency.first;
  return Yield*(DCSCorrection/SingleTagYield);
}

const uncertainties::udouble& FlavourTagYield::GetDTFlavourYield(int Bin) const {
  const std::size_t NumberBins = m_DTFlavourYields.size()/2;
  if(Bin > 0) {
    return m_DTFlavourYields[NumberBins + Bin - 1];
  } else if(Bin < 0) {
    return m_DTFlavourYields[Bin + NumberBins];
  } else {
    throw std::invalid_argument("Bin number must be non-zero");
  }
}

double FlavourTagYield::CalculateDCSCorrection(int Bin,
					       double Ki,
					       double Kbari,
					       const std::vector<double> &ci,
					       const std::vector<double> &si) const {
  if(Bin == 0) {
    throw std::invalid_argument("Bin number must be non-zero");
  }
  const int Sign = Bin > 0 ? +1 : -1;
  const double DCS_term = m_rD.first*m_rD.first*Kbari/Ki;
  const int BinIndex = TMath::Abs(Bin) - 1;
  const double Eff_si = si[BinIndex]*Sign;
  const double Interference_term = -2.0*m_rD.first*m_R.first*(m_CosDeltaD*ci[BinIndex]
                                                            - m_SinDeltaD*Eff_si);
  const double SqrtKK = TMath::Sqrt(Kbari/Ki);
  return 1.0/(1.0 + DCS_term + SqrtKK*Interference_term);
}

double FlavourTagYield::CalculateDCSCorrection(int Bin,
					       const std::vector<double> &ci,
					       const std::vector<double> &si) const {
  double Ki = uncertainties::nom(GetDTFlavourYield(Bin));
  double Kbari = uncertainties::nom(GetDTFlavourYield(-Bin));
  double Ki_temp = Ki;
  double Kbari_temp = Kbari;
  for(std::size_t i = 0; i < 4; i++) {
    Ki_temp = CalculateDCSCorrection(Bin, Ki, Kbari, ci, si);
    Kbari_temp = CalculateDCSCorrection(-Bin, Kbari, Ki, ci, si);
    Ki = Ki_temp;
    Kbari = Kbari_temp;
  }
  return Ki;
}

std::vector<uncertainties::udouble>
FlavourTagYield::LoadFlavourYields(const std::string &TagMode,
				   const Settings &settings) const {
  const std::size_t NumberBins = settings["BinningScheme"].getI("NumberBins");
  std::vector<double> FlavourYields(2*NumberBins);
  std::size_t Index = 0;
  // Loop over all bins
  for(int Bin = -NumberBins; Bin <= static_cast<int>(NumberBins); Bin++) {
    if(Bin == 0) {
      continue;
    }
    const std::string NamePrefix = "DoubleTag_Flavour_KKpipi_vs_" + TagMode +
                                 + "_SignalBin" + (Bin > 0 ? "P" : "M");
    const std::string BinName = std::to_string(TMath::Abs(Bin));
    const std::string VarName = NamePrefix + BinName + "_TagBin0_SignalYield";
    // Get the raw yields
    FlavourYields[Index] = settings[TagMode + "_DT_Yield"].getD(VarName);
    Index++;
  }
  // Load the covariance matrix
  const TMatrixT<double> CovMatrix = LoadCovarianceMatrix(TagMode, settings);
  // Flatten the covariance matrix (for the uncertainties::udouble type)
  std::vector<double> CovMatrix_flat;
  for(std::size_t i = 0; i < 2*NumberBins; i++) {
    for(std::size_t j = 0; j < 2*NumberBins; j++) {
      CovMatrix_flat.push_back(CovMatrix(i, j));
    }
  }
  // Load the raw yields together with their covariance matrix
  std::vector<uncertainties::udouble> FlavourYieldsWithCov =
    uncertainties::ureals<std::vector<uncertainties::udouble>>(FlavourYields, CovMatrix_flat);
  // Load efficiency matrix
  TMatrixT<double> EffMatrix = LoadEfficiencyMatrix(TagMode, settings);
  EffMatrix.Invert();
  // Multiply raw yields with efficiency matrix to obtain efficiency corrected yields
  std::vector<uncertainties::udouble> FlavourYieldsEffCorrected(2*NumberBins);
  for(std::size_t i = 0; i < 2*NumberBins; i++) {
    for(std::size_t j = 0; j < 2*NumberBins; j++) {
      FlavourYieldsEffCorrected[i] += EffMatrix(i, j)*FlavourYieldsWithCov[j];
    }
  }
  return FlavourYieldsEffCorrected;
}

TMatrixT<double>
FlavourTagYield::LoadCovarianceMatrix(const std::string &TagMode,
				      const Settings &settings) const {
  const std::string Filename = settings.get(TagMode + "_RawYields_CorrelationMatrix");
  TFile File(Filename.c_str(), "READ");
  TMatrixTSym<double> *CovMatrix = nullptr;
  File.GetObject("CovarianceMatrix", CovMatrix);
  return TMatrixT<double>(*CovMatrix);
}

TMatrixT<double>
FlavourTagYield::LoadEfficiencyMatrix(const std::string &TagMode,
				      const Settings &settings) const {
  std::string Filename = settings.get(TagMode + "_EfficiencyMatrix");
  TFile File(Filename.c_str(), "READ");
  TMatrixT<double> *EffMatrix = nullptr;
  File.GetObject("EffMatrix", EffMatrix);
  return *EffMatrix;
}

std::pair<double, double>
FlavourTagYield::LoadCharmParameter(const std::string &TagMode,
				    const Settings &settings,
				    const std::string &ParameterName) const {
  const std::string Name = TagMode + "_DCS_Parameters";
  const double Value = settings[Name].getD(ParameterName);
  const double Error = settings[Name].getD(ParameterName + "_err");
  return std::make_pair(Value, Error);
}

TMatrixT<double>
FlavourTagYield::LoadCharmCovMatrix(const std::string &TagMode,
				    const Settings &settings) const {
  const std::string Name = TagMode + "_DCS_Parameters";
  std::string ListNumbers = settings[Name].get("CorrelationMatrix");
  std::replace(ListNumbers.begin(), ListNumbers.end(), ',', ' ');
  std::stringstream ss(ListNumbers);
  const std::vector<double> Uncertainties{m_rD.second, m_R.second, m_DeltaD.second};
  TMatrixT<double> CovMatrix(3, 3);
  for(std::size_t i = 0; i < 3; i++) {
    for(std::size_t j = 0; j < 3; j++) {
      double Correlation;
      ss >> Correlation;
      CovMatrix(i, j) = Correlation*Uncertainties[i]*Uncertainties[j];
    }
  }
  return CovMatrix;
}

std::pair<double, double>
FlavourTagYield::GetSTYield(const std::string &Tag,
			    const Settings &settings) const {
  const std::string SettingsName(Tag + "_ST_Yield");
  const std::string YieldName(Tag + "_SingleTag_Yield");
  const double Yield = settings[SettingsName].getD(YieldName);
  const double Yield_err = settings[SettingsName].getD(YieldName + "_err");
  return std::make_pair(Yield, Yield_err);
}

std::pair<double, double>
FlavourTagYield::GetSTEfficiency(const std::string &Tag,
				 const Settings &settings) const {
  const std::string EffName = Tag + "_SingleTagEfficiency";
  const double Efficiency = settings["ST_Efficiency"].getD(EffName);
  const double Efficiency_err = settings["ST_Efficiency"].getD(EffName + "_err");
  return std::make_pair(Efficiency, Efficiency_err);
}

void FlavourTagYield::PrintDCS(const std::vector<double> &ci,
			  const std::vector<double> &si) const {
  std::cout << std::left << std::setw(10) << "Bin";
  std::cout << std::left << std::setw(10) << m_TagMode + " DCS corrections";
  std::cout << "\n";
  const std::size_t Size = ci.size();
  for(int Bin = -Size; Bin <= static_cast<int>(Size); Bin++) {
    if(Bin == 0) {
      continue;
    }
    std::cout << std::left << std::setw(10) << Bin;
    std::cout << std::left << std::setw(10) << CalculateDCSCorrection(Bin, ci, si);
    std::cout << "\n";
  }
}
