// Martin Duy Tat 3rd October 2022
/**
 * RawBinnedDTYields is an abstract class that stores the binned double tag yields
 * and the corresponding correlation matrix
 */

#ifndef RAWBINNEDDTYIELDS
#define RAWBINNEDDTYIELDS

#include<vector>
#include<unordered_map>
#include<memory>
#include"TMatrixTSym.h"

struct AsymmetricUncertainty {
  /**
   * The central value
   */
  double Value;
  /**
   * The positive uncertainty
   */
  double PlusUncertainty;
  /**
   * The negative uncertainty
   */
  double NegativeUncertainty;
};

class RawBinnedDTYields {
 public:
  /**
   * We don't need the copy constructor
   */
  RawBinnedDTYields(const RawBinnedDTYields &DTYields) = delete;
  /**
   * We don't need the move constructor
   */
  RawBinnedDTYields(RawBinnedDTYields &&DTYields) = delete;
 protected:
  /**
   * Constructor that saves the yields, asymmetric uncertainties and correlation matrix
   * The constructor is protected because it shouldn't exist on its own
   * @param Yields The bin yields with asymmetric uncertainties
   */
  RawBinnedDTYields(const std::vector<AsymmetricUncertainty> &Yields,
		    const TMatrixTSym<double> &CorrelationMatrix);
 public:
  /**
   * Get the double tag yields
   */
  const std::vector<AsymmetricUncertainty>& GetDoubleTagYields() const;
  /**
   * Get the correlation matrix
   */
  const TMatrixTSym<double>& GetCorrelationMatrix() const;
  /**
   * Generate toy yields with a uniform distribution in a 5 sigma window
   */
  std::vector<AsymmetricUncertainty> GetToyYields() const;
 private:
  /**
   * The binned double tag yields with asymmetric uncertainties
   */
  const std::vector<AsymmetricUncertainty> m_Yields;
  /**
   * The correlation matrix
   */
  const TMatrixTSym<double> m_CorrelationMatrix;
};

#endif
