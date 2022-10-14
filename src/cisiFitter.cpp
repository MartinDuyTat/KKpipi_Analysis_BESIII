// Martin Duy Tat 13th October 2022

#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"Settings.h"
#include"cisiLikelihood.h"
#include"cisiFitter.h"

cisiFitter::cisiFitter(const Settings &settings):
  m_cisiLikelihood(settings),
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")) {
}

void cisiFitter::Minimise() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(4);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    std::vector<double> ci, si;
    for(int i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i]);
      si.push_back(x[i + m_NumberBins]);
    }
    return m_cisiLikelihood.CalculateLogLikelihood(ci, si);
  };
  ROOT::Math::Functor fcn(LikelihoodFunction, 2*m_NumberBins);
  Minimiser.SetFunction(fcn);
  for(int i = 0; i < m_NumberBins; i++) {
    Minimiser.SetVariable(i, "c" + std::to_string(i + 1), 1.0, 1.0);
    Minimiser.SetVariableLimits(i, -2.0, 2.0);
  }
  for(int i = 0; i < m_NumberBins; i++) {
    Minimiser.SetVariable(i + m_NumberBins, "s" + std::to_string(i + 1), 0.0, 1.0);
    Minimiser.SetVariableLimits(i + m_NumberBins, -2.0, 2.0);
    Minimiser.FixVariable(i + m_NumberBins);
  }
  Minimiser.Minimize();
  const double *X = Minimiser.X();
  std::vector<double> ci;
  std::vector<double> si;
  for(int i = 0; i < m_NumberBins; i++) {
    ci.push_back(X[i]);
    si.push_back(X[i + m_NumberBins]);
  }
  m_cisiLikelihood.PrintComparison(ci, si);
}
