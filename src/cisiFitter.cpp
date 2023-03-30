// Martin Duy Tat 13th October 2022

#include"TFile.h"
#include"TTree.h"
#include"TRandom.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"
#include"TLatex.h"
#include"TLine.h"
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
    for(std::size_t i = 0; i < m_NumberBins; i++) {
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
  for(std::size_t i = 0; i < m_NumberBins; i++) {
    ci.push_back(X[i + 1]);
    si.push_back(X[i + m_NumberBins + 1]);
  }
  m_cisiLikelihood.PrintComparison(X[0], ci, si);
  m_cisiLikelihood.PrintKi(ci, si);
  Plot_cisi(Minimiser, ci, si);
}

void cisiFitter::RunToys() const {
  gRandom->SetSeed(42);
  const double Generator_BF_KKpipi = 0.00247;
  double BF_KKpipi_value, BF_KKpipi_err, BF_KKpipi_pull;
  int Status, CovStatus;
  double LL;
  const std::vector<double> Generator_ci = GetGeneratorcisi("c");
  const std::vector<double> Generator_si = GetGeneratorcisi("s");
  std::vector<double> ci_value(m_NumberBins), si_value(m_NumberBins);
  std::vector<double> ci_err(m_NumberBins), si_err(m_NumberBins);
  std::vector<double> ci_pull(m_NumberBins), si_pull(m_NumberBins);
  TFile File(m_Settings.get("ManyToysOutputFilename").c_str(), "RECREATE");
  TTree Tree("cisiTree", "");
  Tree.Branch("Status", &Status);
  Tree.Branch("CovStatus", &CovStatus);
  Tree.Branch("LL", &LL);
  Tree.Branch("BF_KKpipi_value", &BF_KKpipi_value);
  Tree.Branch("BF_KKpipi_err", &BF_KKpipi_err);
  Tree.Branch("BF_KKpipi_pull", &BF_KKpipi_pull);
  for(std::size_t i = 0; i < m_NumberBins; i++) {
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
    for(std::size_t i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i + 1]);
      si.push_back(x[i + m_NumberBins + 1]);
    }
    return m_cisiLikelihood.CalculateToyLogLikelihood(BF_KKpipi, ci, si);
  };
  std::size_t StatsMultiplier = m_Settings.getI("StatsMultiplier");
  std::size_t NumberOfToysToSkip = m_Settings.getI("NumberOfToysToSkip");
  for(std::size_t i = 0; i < NumberOfToysToSkip; i++) {
    std::cout << "Skipped toy " << i << "\n";
    cisiLikelihoodRef.GenerateToy(Generator_BF_KKpipi,
				  Generator_ci,
				  Generator_si,
				  StatsMultiplier);
  }
  std::size_t NumberToys = m_Settings.getI("NumberToys");
  for(std::size_t i = 0; i < NumberToys; i++) {
    std::cout << "Toy " << i << "\n";
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
    LL = Minimiser.MinValue();
    const double *X = Minimiser.X();
    const double *E = Minimiser.Errors();
    BF_KKpipi_value = X[0];
    BF_KKpipi_err = E[0];
    BF_KKpipi_pull = (BF_KKpipi_value - Generator_BF_KKpipi)/BF_KKpipi_err;
    for(std::size_t i = 0; i < m_NumberBins; i++) {
      ci_value[i] = X[i + 1];
      si_value[i] = X[i + m_NumberBins + 1];
      ci_err[i] = E[i + 1];
      si_err[i] = E[i + m_NumberBins + 1];
      double ErrorLow, ErrorHigh;
      Minimiser.GetMinosError(i + 1, ErrorLow, ErrorHigh, 0);
      ci_pull[i] = ci_value[i] - Generator_ci[i];
      //ci_pull[i] /= ci_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      ci_pull[i] /= ci_err[i];
      si_pull[i] = si_value[i] - Generator_si[i];
      Minimiser.GetMinosError(i + 1 + m_NumberBins, ErrorLow, ErrorHigh, 0);
      //si_pull[i] /= si_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      si_pull[i] /= si_err[i];
    }
    Tree.Fill();
  }
  Tree.Write();
  File.Close();
}

void cisiFitter::SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const {
  Minimiser.SetVariable(0, "BF_KKpipi", 0.00247, 0.1);
  Minimiser.SetVariableLimits(0, 0.001, 0.010);
  const double cMin = m_Settings.getD("c_min");
  const double cMax = m_Settings.getD("c_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "c" + std::to_string(i);
    Minimiser.SetVariable(i, Name, 1.0, 1.0);
    Minimiser.SetVariableLimits(i, cMin, cMax);
  }
  const double sMin = m_Settings.getD("s_min");
  const double sMax = m_Settings.getD("s_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "s" + std::to_string(i);
    Minimiser.SetVariable(i + m_NumberBins, Name, 0.0, 1.0);
    Minimiser.SetVariableLimits(i + m_NumberBins, sMin, sMax);
    if(m_Settings.getB("Fix_si")) {
      Minimiser.FixVariable(i + m_NumberBins);
    }
  }
}

std::vector<double> cisiFitter::GetGeneratorcisi(const std::string &c_or_s) const {
  std::vector<double> cisi(m_NumberBins);
  const std::string Prefix = "Generator_" + c_or_s;
  for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
    cisi[Bin - 1] = m_Settings.getD(Prefix + std::to_string(Bin));
  }
  return cisi;
}

void cisiFitter::Plot_cisi(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
			   const std::vector<double> &ci,
			   const std::vector<double> &si) const {
  // Make canvas
  TCanvas c("cisi_c", "", 900, 900);
  // Draw circle
  const std::size_t Points = 101;
  std::vector<double> circle_x(Points), circle_y(Points);
  for(std::size_t i = 0; i < 101; i++) {
    circle_x[i] = TMath::Cos(2.0*TMath::Pi()*i/100.0);
    circle_y[i] = TMath::Sin(2.0*TMath::Pi()*i/100.0);
  }
  TGraph Circle(Points, circle_x.data(), circle_y.data());
  Circle.SetLineWidth(3);
  Circle.SetTitle(";#it{c_{i}};#it{s_{i}}");
  Circle.Draw("AL");
  const double Boundary = m_Settings.getD("PlotBoundary");
  Circle.GetXaxis()->SetLimits(-Boundary, Boundary);
  Circle.GetYaxis()->SetRangeUser(-Boundary, Boundary);
  Circle.GetXaxis()->SetNdivisions(505);
  Circle.GetYaxis()->SetNdivisions(505);
  // Draw prediction of F+
  const double c_fromFPlus = 2*0.73 - 1;
  const double c_fromFPlus_err = 2*0.04;
  TLine LeftLine(c_fromFPlus - c_fromFPlus_err, -Boundary,
		 c_fromFPlus - c_fromFPlus_err, Boundary);
  TLine RightLine(c_fromFPlus + c_fromFPlus_err, -Boundary,
		  c_fromFPlus + c_fromFPlus_err, Boundary);
  LeftLine.SetLineWidth(3);
  RightLine.SetLineWidth(3);
  LeftLine.SetLineStyle(kDashed);
  RightLine.SetLineStyle(kDashed);
  LeftLine.Draw("L SAME");
  RightLine.Draw("L SAME");
  // Draw fit results
  TGraph Results(m_NumberBins, ci.data(), si.data());
  Results.SetMarkerStyle(8);
  std::vector<TLatex> PointLabels;
  for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
    PointLabels.push_back(TLatex(ci[Bin - 1], si[Bin - 1] + 0.03, std::to_string(Bin).c_str()));
    PointLabels.back().DrawClone("SAME");
  }
  Results.Draw("P SAME");
  // Create and draw contours
  Minimiser.SetPrintLevel(-1);
  std::vector<TGraph> Contours;
  std::vector<double> ErrorDefs{1.0, 4.0};
  for(double ErrorDef : ErrorDefs) {
    Minimiser.SetErrorDef(ErrorDef);
    for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
      std::vector<double> x(Points), y(Points);
      unsigned int nPoints = Points - 1;
      Minimiser.Contour(Bin, Bin + m_NumberBins, nPoints, x.data(), y.data());
      x.back() = x[0];
      y.back() = y[0];
      Contours.push_back(TGraph(Points, x.data(), y.data()));
      Contours.back().SetLineWidth(3);
    }
  }
  for(auto &Contour : Contours) {
    Contour.Draw("L SAME");
  }
  c.SaveAs("Contours_cisi.pdf");
}
