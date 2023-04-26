// Martin Duy Tat 3rd October 2022
/**
 * RawBinnedFlavourTagYields inhertis from RawBinnedDTYields
 * It initialises binned yields of a flavour tag
 */

#ifndef RAWBINNEDFLAVOURTAGYIELDS
#define RAWBINNEDFLAVOURTAGYIELDS

#include<string>
#include<vector>
#include"RawBinnedDTYields.h"
#include"Settings.h"

class RawBinnedFlavourTagYields: public RawBinnedDTYields {
 public:
  /**
   * Constructor that parses the fitted double tag yields
   * @param Tag The name of the tag
   * @param settings The settings file
   * @param ToyNumber For toys, this number labels which toy, otherwise set to -1
   */
  RawBinnedFlavourTagYields(const std::string &Tag,
			    const Settings &settings,
			    int ToyNumber = -1);
 private:
  /**
   * Helper function that parses the binned double tag yields from a file
   * @param Tag The tag mode
   * @param settings The settings file
   * @param ToyNumber For toys, this number labels which toy, otherwise set to -1
   */
  std::vector<AsymmetricUncertainty>
  ParseYields(const std::string &Tag,
	      const Settings &settings,
	      int ToyNumber) const;
};

#endif
