// Martin Duy Tat 3rd October 2022

#include<vector>
#include"TMatrixTSym.h"
#include"RawBinnedDTYields.h"

RawBinnedDTYields::RawBinnedDTYields(const std::vector<AsymmetricUncertainty> &Yields,
				     const TMatrixTSym<double> &CorrelationMatrix):
  m_Yields(Yields),
  m_CorrelationMatrix(CorrelationMatrix) {
}

const std::vector<AsymmetricUncertainty>&
RawBinnedDTYields::GetDoubleTagYields() const {
  return m_Yields;
}

const TMatrixTSym<double>& RawBinnedDTYields::GetCorrelationMatrix() const {
  return m_CorrelationMatrix;
}
