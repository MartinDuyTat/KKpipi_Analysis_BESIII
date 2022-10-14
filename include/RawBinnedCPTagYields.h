// Martin Duy Tat 3rd October 2022
/**
 * RawBinnedCPTagYields inhertis from RawBinnedDTYields
 * It initialises binned yields of a CP tag
 */

#include<string>

#ifndef RAWBINNEDCPTAGYIELDS
#define RAWBINNEDCPTAGYIELDS

#include<string>
#include<vector>
#include"RawBinnedDTYields.h"
#include"Settings.h"

class RawBinnedCPTagYields: public RawBinnedDTYields {
 public:
  /**
   * Constructor that parses the fitted double tag yields
   * @param Tag The name of the tag
   * @param settings The settings file
   */
  RawBinnedCPTagYields(const std::string &Tag, const Settings &settings);
 private:
  /**
   * Helper function that parses the binned double tag yields from a file
   * @param Tag The tag mode
   * @param settings The settings file
   */
  std::vector<AsymmetricUncertainty>
  ParseYields(const std::string &Tag, const Settings &settings) const;
  /**
   * Helper function to load the correlation matrix
   */
  TMatrixTSym<double> LoadCorrelationMatrix(const std::string &Tag,
					    const Settings &settings) const;
};

#endif
