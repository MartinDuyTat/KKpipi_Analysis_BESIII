// Martin Duy Tat 31st March 2021
/**
 * Utilities is a namespace with useful utility functions
 */

#include<string>
#include"TChain.h"

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
}

#endif
