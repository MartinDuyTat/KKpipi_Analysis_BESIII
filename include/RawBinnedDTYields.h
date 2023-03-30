// Martin Duy Tat 3rd October 2022
/**
 * RawBinnedDTYields is an abstract class that stores the binned double tag yields
 * and the corresponding correlation matrix
 */

#ifndef RAWBINNEDDTYIELDS
#define RAWBINNEDDTYIELDS

#include<vector>
#include<memory>
#include<utility>
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
  /**
   * The symmetric uncertainty (Hessian)
   */
  double SymmetricUncertainty;
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
   * Generate toy yields, drawn from a Gaussian distribution
   * @param PredictedYields The predicted yields based on input parameters
   * @param SymmetricUncertainties Set to true for symmetric uncertainties
   * @return The toy yields and the probability of the underlying Gaussian distribution
   */
  std::pair<std::vector<AsymmetricUncertainty>, double>
  GetToyYields(const std::vector<double> &PredictedYields,
	       bool SymmetricUncertainties) const;
 private:
  /**
   * The binned double tag yields with asymmetric uncertainties
   */
  const std::vector<AsymmetricUncertainty> m_Yields;
  /**
   * The correlation matrix
   */
  const TMatrixTSym<double> m_CorrelationMatrix;
  /**
   * Helper function that generates a single yield from asymmetric uncertainties
   * @param PredictedYield The predicted yield given ci and si
   * @param DataYield The yield and uncertainties measured in data
   * @param SymmetricUncertainties Set to true to use the symmetric uncertaintnies
   */
  std::pair<double, double> GenerateYield(double PredictedYield,
					  const AsymmetricUncertainty &DataYield,
					  bool SymmetricUncertainties) const;
};

#endif
