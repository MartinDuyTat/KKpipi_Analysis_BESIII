// Martin Duy Tat 3rd October 2022

#include<vector>
#include<memory>
#include"TMatrixTSym.h"
#include"TRandom.h"
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

std::vector<AsymmetricUncertainty> RawBinnedDTYields::GetToyYields() const {
  auto ToyYields = m_Yields;
  for(std::size_t i = 0; i < m_Yields.size(); i++) {
    const double LowerLimit = m_Yields[i].Value
                            - 5.0*m_Yields[i].NegativeUncertainty;
    const double UpperLimit = m_Yields[i].Value
                            + 5.0*m_Yields[i].PlusUncertainty;
    ToyYields[i].Value = gRandom->Uniform(LowerLimit, UpperLimit);
  }
  return ToyYields;
}
