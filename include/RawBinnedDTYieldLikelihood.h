// Martin Duy Tat 12th May 2023
/**
 * RawBinnedDTYieldLikelihood contains the full double tag yield likelihood
 */

#ifndef RAWBINNEDDTYIELDLIKELIHOOD
#define RAWBINNEDDTYIELDLIKELIHOOD

#include<vector>
#include<string>
#include<memory>
#include"RooArgSet.h"
#include"RooAbsReal.h"
#include"Settings.h"

class RawBinnedDTYieldLikelihood {
 public:
  /**
   * Constructor that loads the likelihood from a file
   * @param Tag The name of the tag
   * @param settings The settings file
   * @param TagCategory CP, Flavour or SCMB
   * @param ToyNumber For toys, this number labels which toy, otherwise set to -1
   */
  RawBinnedDTYieldLikelihood(const std::string &Tag,
			     const Settings &settings,
			     const std::string &TagCategory,
			     int ToyNumber = -1);
  /**
   * Constructor that loads the Feldman Cousins likelihood from a file
   * @param Tag The name of the tag
   * @param settings The settings file
   * @param TagCategory CP, Flavour or SCMB
   * @param ToyName The name of the toy
   * @param ToyNumber Feldman Cousins toy number
   */
  RawBinnedDTYieldLikelihood(const std::string &Tag,
			     const Settings &settings,
			     const std::string &TagCategory,
			     const std::string &ToyName,
			     int ToyNumber);
  /**
   * We don't need copy constructor
   */
  RawBinnedDTYieldLikelihood(const RawBinnedDTYieldLikelihood &LL) = delete;
  /**
   * We own the pointer to the likelihood, need to delete
   */
  ~RawBinnedDTYieldLikelihood() = default;
  /**
   * Evaluate the log likelihood
   */
  double GetLogLikelihood(const std::vector<double> &PredictedBinYields) const;
 private:
  /**
   * Tag mode
   */
  const std::string m_TagMode;
  /**
   * The tag category (CP, Flavour or SCMB)
   */
  const std::string m_TagCategory;
  /**
   * The full likelihood of the signal yields
   */
  std::unique_ptr<RooAbsReal> m_FullLikelihood;
  /**
   * The fit variables
   */
  std::unique_ptr<RooArgSet> m_Variables;
  /**
   * The ordering of the yields
   */
  std::vector<std::string> m_Order;
  /**
   * Helper function that loads the full likelihood from a file
   */
  RooAbsReal* GetFullLikelihood(const std::string &Tag,
				const Settings &settings,
				int ToyNumber) const;
  /**
   * Helper function to set up the yield ordering
   */
  std::vector<std::string> GetYieldOrder(int NumberBins) const;
};

#endif
