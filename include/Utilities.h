// Martin Duy Tat 31st March 2021
/**
 * Utilities is a namespace with useful utility functions
 */

#include<string>
#include"TChain.h"
#include"TCut.h"

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
   * @param CutType "DeltaECuts" for standard initial cuts plus \f$\Delta E\f$ cuts, "NoDeltaECuts" for standard initial cuts only, "TruthMatchingCuts" for truth matching cuts
   * @param TagMode "KKpipi", "Kpi", etc
   * @param TagType "ST" for single tag and "DT" for double tag
   * @param DataMC "Data" or "MC"
   */
  TCut LoadCuts(const std::string &CutType, const std::string &TagMode, const std::string &TagType, const std::string &DataMC);
}

#endif
