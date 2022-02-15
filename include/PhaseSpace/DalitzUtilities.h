// Martin Duy Tat 14th February 2022
/**
 * DalitzUtilities is a namespace containing useful functions to analyse the Dalitz phase space
 */

#ifndef DALITZUTILITIES
#define DALITZUTILITIES

#include"TH2F.h"

namespace DalitzUtilities {
  /**
   * Read bin number from lookup table
   */
  int LookUpBinNumber(double M2Plus, double M2Minus, const TH2F *BinningScheme);
  /**
   * Map a Dalitz point outside of phase space back inside the boundary
   */
  int GetMappedK0hhBin(double M2Plus, double M2Minus, const TH2F *BinningScheme);
}

#endif
