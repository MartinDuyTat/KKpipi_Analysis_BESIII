// Martin Duy Tat 12th May 2023
/**
 * RawBinnedDTYieldLikelihood contains the full double tag yield likelihood
 */

#ifndef RAWBINNEDDTYIELDLIKELIHOOD
#define RAWBINNEDDTYIELDLIKELIHOOD

#include<vector>
#include<string>
#include"RooNLLVar.h"
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
   * We don't need copy constructor
   */
  RawBinnedDTYieldLikelihood(const RawBinnedDTYieldLikelihood &LL) = delete;
  /**
   * We own the pointer to the likelihood, need to delete
   */
  ~RawBinnedDTYieldLikelihood();
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
  RooNLLVar *m_FullLikelihood;
  /**
   * The constant offset at the minimum
   */
  double m_Offset;
  /**
   * The ordering of the yields
   */
  std::vector<std::string> m_Order;
  /**
   * Helper function that loads the full likelihood from a file
   */
  RooNLLVar* GetFullLikelihood(const std::string &Tag,
			       const Settings &settings,
			       int ToyNumber) const;
  /**
   * Helper function to set up the yield ordering
   */
  std::vector<std::string> GetYieldOrder(int NumberBins) const;
};

#endif
