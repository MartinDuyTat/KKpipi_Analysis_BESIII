// Martin Duy Tat 31st March 2021
/**
 * Utilities is a namespace with useful utility functions
 */

#include<string>
#include<vector>
#include<memory>
#include<numeric>
#include"TChain.h"
#include"TTree.h"
#include"TCut.h"
#include"RooRealVar.h"
#include"Settings.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"
#include"PhaseSpace/KKpipi_vs_Flavour_PhaseSpace.h"

#ifndef UTILITIES
#define UTILITIES

namespace Utilities {
  /**
   * LoadChain is a function that takes in a TChain pointer and loads the ROOT files with TTree objects
   * ROOT files should have a filename ending with a number, starting at 0, before .root
   * @param Chain TChain pointer where TTree objects are loaded
   * @param NumberFiles Number of files
   * @param Filename Filename of ROOT files, without the number and .root
   * @param TreeName Name of the TTree we want to load, if this is not given it is assumed the name of the TChain has already been set
   */
  void LoadChain(TChain *Chain, int NumberFiles, const std::string &Filename, const std::string &TreeName = std::string());
  /**
   * This is an overload which will load a TChain by reading the filenames from a text file
   * @param Chain TChain pointer where TTree objects are loaded
   * @param Filename Filename of text file containing all the individual ROOT files
   * @param TreeName Name of the TTree we want to load, if this is not given it is assumed the name of the TChain has already been set
   */
  void LoadChain(TChain *Chain, const std::string &Filename, const std::string &TreeName = std::string());
  /**
   * This is a helper function for loading the correct cuts
   * @param SignalMode "KKpipi", "KSKK_to_KKpipi", etc
   * @param TagMode "KKpipi", "Kpi", "Kpi_to_KK", etc
   * @param TagType "ST" for single tag and "DT" for double tag
   * @param IncludeDeltaECuts Set to true to apply \f$Delta E\f$ cuts
   * @param DataMC "Data" or "MC"
   * @param TruthMatch Set to true to truth match
   * @param KKpipiPartReco Set to true if KKpipi signal mode is partially reconstructed
   */
  TCut LoadCuts(const std::string &SignalMode, const std::string &TagMode, const std::string &TagType, const std::string &DataMC, bool IncludeDeltaECuts, bool TruthMatch = false, bool KKpipiPartReco = false);
  /**
   * Parse arguments and set up a Settings object
   * Copied from GGSZ code repository
   * Pass inputs to application
   */
  Settings parse_args(int argc, char** argv);
  /**
   * Helper function for Settings class
   * Copied from GGSZ code repository
   */
  void replace_env_variables(std::string& text);
  /**
   * Load a parameter and initialize the RooRealVar object
   * Copied from GGSZ code repository
   * @param settings The settings containing the parameters
   * @param name Name of the RooRealVar variable
   */
  RooRealVar* load_param(const Settings &settings, const std::string &name);
  /**
   * Convert a comma separated list into a vector
   * @param List List of comma separated elements
   */
  std::vector<std::string> ConvertStringToVector(std::string List);
  /**
   * Replace part of string with another string
   * @param String The original string
   * @param From The string we want to replace
   * @param To The string we're replacing with
   */
  std::string ReplaceString(const std::string &String, const std::string &From, const std::string &To);
  /**
   * Get the correct phase space binning object
   * @param settings The fit settings
   * @param Tree The tree containing double tag events to be binned
   */
  std::unique_ptr<KKpipi_PhaseSpace> GetPhaseSpaceBinning(const Settings &settings, TTree *Tree);
  /**
   * Determine the tag type, which can be "Flavour", "CP", "SCMB"
   */
  std::string GetTagType(const std::string &Mode);
  /**
   * Calculate covariance of two vectors
   */
  template<typename T>
  T Covariance(const std::vector<T> &x, const std::vector<T> &y) {
    T Mean_x = TMath::Mean(x.begin(), x.end());
    T Mean_y = TMath::Mean(y.begin(), y.end());
    T Total2 = std::inner_product(x.begin(), x.end(), y.begin(), static_cast<T>(0),
				  std::plus<>(), [=](T a, T b) { return (a - Mean_x)*(b - Mean_y); });
    return Total2/static_cast<T>(x.size() - 1);
  }
  /**
   * Sum the weights of all events that pass the cut
   * Tree The tree with events
   * Cut The cut applied
   */
  double SumWeights(TTree *Tree, const std::string &WeightName, const std::string &Cut = "");
  /**
   * Get the correct ROOT LaTeX name for the tag mode
   */
  std::string GetTagNameLaTeX(const std::string &Tag);
}

#endif
