// Martin Duy Tat 21st April 2022
/**
 * Cholesky smearing is a class for smearing of correlated parameters
 * The smearing is generated from the covariance matrix by Choleksy decomposing the covariance matrix
 * The smearings are saved inside the class and can be retrieved with a getter
 */

#ifndef CHOLESKYSMEARING
#define CHOLESKYSMEARING

#include"TMatrixT.h"

class CholeskySmearing {
  public:
    /**
     * Constructor that initializes the smearing matrices
     * @param CovMatrix Covariance matrix of parameters we want to smear
     */
    CholeskySmearing(const TMatrixT<double> &CovMatrix);
    /**
     * Generate new smearing parameters
     */
    void Smear();
    /**
     * Get smeared parameter
     * @param i Index labelling the parameter
     */
    double GetSmearing(int i) const;
    /**
     * Get the whole vector of smeared parameters
     */
    TMatrixT<double> GetSmearings() const;
  private:
    /**
     * The Cholesky decomposition of the covariance matrix
     */
    const TMatrixT<double> m_CholeskyMatrix;
    /**
     * The smearing of the parameters
     */
    TMatrixT<double> m_Smearings;
    /**
     * Helper function for obtaining the Cholesky decomposition of the covariance matrix
     */
    TMatrixT<double> GetCholeskyDecomposition(const TMatrixT<double> &CovMatrix) const;
};

#endif
