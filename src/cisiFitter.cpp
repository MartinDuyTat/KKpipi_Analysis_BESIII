// Martin Duy Tat 13th October 2022

#include"TFile.h"
#include"TTree.h"
#include"TRandom.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"Settings.h"
#include"cisiLikelihood.h"
#include"cisiFitter.h"

cisiFitter::cisiFitter(const Settings &settings):
  m_cisiLikelihood(settings),
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")),
  m_Settings(settings) {
}

void cisiFitter::Minimise() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(4);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const double BF_KKpipi = x[0];
    std::vector<double> ci, si;
    for(int i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i + 1]);
      si.push_back(x[i + m_NumberBins + 1]);
    }
    return m_cisiLikelihood.CalculateLogLikelihood(BF_KKpipi, ci, si);
  };
  ROOT::Math::Functor fcn(LikelihoodFunction, 2*m_NumberBins + 1);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  Minimiser.Minimize();
  const double *X = Minimiser.X();
  std::vector<double> ci;
  std::vector<double> si;
  for(int i = 0; i < m_NumberBins; i++) {
    ci.push_back(X[i + 1]);
    si.push_back(X[i + m_NumberBins + 1]);
  }
  m_cisiLikelihood.PrintComparison(X[0], ci, si);
}

void cisiFitter::RunToys() const {
  gRandom->SetSeed(42);
  const double Generator_BF_KKpipi = 0.00247;
  double BF_KKpipi_value, BF_KKpipi_err, BF_KKpipi_pull;
  int Status, CovStatus;
  const std::vector<double> Generator_ci{0.50205, 0.588241}, Generator_si{-0.41232, 0.394319};
  std::vector<double> ci_value(m_NumberBins), si_value(m_NumberBins);
  std::vector<double> ci_err(m_NumberBins), si_err(m_NumberBins);
  std::vector<double> ci_pull(m_NumberBins), si_pull(m_NumberBins);
  TFile File(m_Settings.get("ManyToysOutputFilename").c_str(), "RECREATE");
  TTree Tree("cisiTree", "");
  Tree.Branch("Status", &Status);
  Tree.Branch("CovStatus", &CovStatus);
  Tree.Branch("BF_KKpipi_value", &BF_KKpipi_value);
  Tree.Branch("BF_KKpipi_err", &BF_KKpipi_err);
  Tree.Branch("BF_KKpipi_pull", &BF_KKpipi_pull);
  for(int i = 0; i < m_NumberBins; i++) {
    Tree.Branch(("c" + std::to_string(i + 1) + "_value").c_str(), &ci_value[i]);
    Tree.Branch(("s" + std::to_string(i + 1) + "_value").c_str(), &si_value[i]);
    Tree.Branch(("c" + std::to_string(i + 1) + "_err").c_str(), &ci_err[i]);
    Tree.Branch(("s" + std::to_string(i + 1) + "_err").c_str(), &si_err[i]);
    Tree.Branch(("c" + std::to_string(i + 1) + "_pull").c_str(), &ci_pull[i]);
    Tree.Branch(("s" + std::to_string(i + 1) + "_pull").c_str(), &si_pull[i]);
  }
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(-1);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const double BF_KKpipi = x[0];
    std::vector<double> ci, si;
    for(int i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i + 1]);
      si.push_back(x[i + m_NumberBins + 1]);
    }
    return m_cisiLikelihood.CalculateToyLogLikelihood(BF_KKpipi, ci, si);
  };
  std::size_t NumberToys = m_Settings.getI("NumberToys");
  std::size_t StatsMultiplier = m_Settings.getI("StatsMultiplier");
  for(std::size_t i = 0; i < NumberToys; i++) {
    cisiLikelihoodRef.GenerateToy(Generator_BF_KKpipi,
				  Generator_ci,
				  Generator_si,
				  StatsMultiplier);
    ROOT::Math::Functor fcn(LikelihoodFunction, 2*m_NumberBins + 1);
    Minimiser.SetFunction(fcn);
    SetupMinimiser(Minimiser);
    Minimiser.Minimize();
    Status = Minimiser.Status();
    CovStatus = Minimiser.CovMatrixStatus();
    const double *X = Minimiser.X();
    const double *E = Minimiser.Errors();
    BF_KKpipi_value = X[0];
    BF_KKpipi_err = E[0];
    BF_KKpipi_pull = (BF_KKpipi_value - Generator_BF_KKpipi)/BF_KKpipi_err;
    for(int i = 0; i < m_NumberBins; i++) {
      ci_value[i] = X[i + 1];
      si_value[i] = X[i + m_NumberBins + 1];
      ci_err[i] = E[i + 1];
      si_err[i] = E[i + m_NumberBins + 1];
      double ErrorLow, ErrorHigh;
      Minimiser.GetMinosError(i + 1, ErrorLow, ErrorHigh, 0);
      ci_pull[i] = ci_value[i] - Generator_ci[i];
      ci_pull[i] /= ci_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      si_pull[i] = si_value[i] - Generator_si[i];
      Minimiser.GetMinosError(i + 1 + m_NumberBins, ErrorLow, ErrorHigh, 0);
      si_pull[i] /= si_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
    }
    Tree.Fill();
  }
  Tree.Write();
  File.Close();
}

void cisiFitter::SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const {
  Minimiser.SetVariable(0, "BF_KKpipi", 0.00247, 0.1);
  Minimiser.SetVariableLimits(0, 0.001, 0.010);
  for(int i = 1; i <= m_NumberBins; i++) {
    Minimiser.SetVariable(i, "c" + std::to_string(i), 1.0, 1.0);
    Minimiser.SetVariableLimits(i, -2.0, 2.0);
  }
  for(int i = 1; i <= m_NumberBins; i++) {
    Minimiser.SetVariable(i + m_NumberBins, "s" + std::to_string(i), 0.0, 1.0);
    Minimiser.SetVariableLimits(i + m_NumberBins, -3.0, 3.0);
    if(m_Settings.getB("Fix_si")) {
      Minimiser.FixVariable(i + m_NumberBins);
    }
  }
}
