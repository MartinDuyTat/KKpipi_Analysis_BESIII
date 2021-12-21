// Martin Duy Tat 31st March 2021
/**
 * Utilities is a namespace with useful utility functions
 */

#include<string>
#include<vector>
#include<memory>
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
   * @param Set to true to truth match
   */
  TCut LoadCuts(const std::string &SignalMode, const std::string &TagMode, const std::string &TagType, const std::string &DataMC, bool IncludeDeltaECuts, bool TruthMatch = false);
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
}

#endif
