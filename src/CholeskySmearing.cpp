// Martin Duy Tat 21st April 2022

#include<stdexcept>
#include"TMatrixT.h"
#include"TRandom.h"
#include"TDecompChol.h"
#include"CholeskySmearing.h"

CholeskySmearing::CholeskySmearing(const TMatrixT<double> &CovMatrix): m_CholeskyMatrix(GetCholeskyDecomposition(CovMatrix)),
								       m_Smearings(CovMatrix.GetNrows(), 1) {
  if(CovMatrix.GetNrows() != CovMatrix.GetNcols()) {
    throw std::range_error("Covariance matrix is not square");
  }
}

void CholeskySmearing::Smear() {
  for(int i = 0; i < m_Smearings.GetNrows(); i++) {
    m_Smearings(i, 0) = gRandom->Gaus(0.0, 1.0);
  }
  m_Smearings = m_CholeskyMatrix*m_Smearings;
}

double CholeskySmearing::GetSmearing(int i) const {
  return m_Smearings(i, 0);
}

TMatrixT<double> CholeskySmearing::GetCholeskyDecomposition(const TMatrixT<double> &CovMatrix) const {
  TDecompChol CholeskyDecomposition(CovMatrix);
  bool Success = CholeskyDecomposition.Decompose();
  if(!Success) {
    throw std::runtime_error("Covariance matrix not positive definite");
  }
  TMatrixT<double> Temp = CholeskyDecomposition.GetU();
  return Temp.T();
}
