// Martin Duy Tat 15th February 2023
/**
 * KiCombiner combines the flavour tag yields of all tags and normalises them
 */

#ifndef KICOMBINER
#define KICOMBINER

#include<vector>
#include"uncertainties/ureal.hpp"
#include"Settings.h"
#include"FlavourTags/FlavourTagYield.h"

class KiCombiner {
 public:
  /**
   * Constructor that initalises all the flavour tags
   */
  KiCombiner(const Settings &settings);
  /**
   * Get the value of Ki, given a value of ci and si
   */
  std::vector<double> GetKi(const std::vector<double> &ci,
			    const std::vector<double> &si) const;
  /**
   * Get the value of Ki, given a value of ci and si, with uncertainties
   */
  std::vector<uncertainties::udouble> GetKiWithUncertainties(
    const std::vector<double> &ci,
    const std::vector<double> &si) const;
  /**
   * Print the values of Ki
   */
  void PrintKi(const std::vector<double> &ci,
	       const std::vector<double> &si) const;
 private:
  /**
   * Number of bins
   */
  const std::size_t m_NumberBins;
  /**
   * Information about each flavour tag
   */
  const std::vector<FlavourTagYield> m_FlavourTags;
  /**
   * Helper function that loads all the flavour tags
   */
  std::vector<FlavourTagYield> LoadFlavourTags(const Settings &settings) const;
};

#endif
