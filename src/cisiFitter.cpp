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
  m_Settings(settings),
  m_FitDeltaKpi(settings.getB("FitDeltaKpi")) {
}

void cisiFitter::Minimise() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(4);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const double BF_KKpipi = x[0];
    std::vector<double> ci, si, Ki, Kbari;
    double DeltaKpi = 0.0;
    for(std::size_t i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i + 0*m_NumberBins + 1]);
      si.push_back(x[i + 1*m_NumberBins + 1]);
      Ki.push_back(x[i + 2*m_NumberBins + 1]);
      Kbari.push_back(x[i + 3*m_NumberBins + 1]);
    }
    if(m_FitDeltaKpi) {
      DeltaKpi = x[4*m_NumberBins + 1];
    }
    return m_cisiLikelihood.CalculateLogLikelihood(BF_KKpipi,
						   ci, si,
						   Ki, Kbari,
						   DeltaKpi);
  };
  std::size_t NumberParameters = 4*m_NumberBins + 1;
  if(m_FitDeltaKpi) {
    NumberParameters++;
  }
  ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  Minimiser.Minimize();
  const double *X = Minimiser.X();
  std::vector<double> ci;
  std::vector<double> si;
  std::vector<double> Ki;
  std::vector<double> Kbari;
  double DeltaKpi = 0.0;
  for(std::size_t i = 0; i < m_NumberBins; i++) {
    ci.push_back(X[i + 0*m_NumberBins + 1]);
    si.push_back(X[i + 1*m_NumberBins + 1]);
    Ki.push_back(X[i + 2*m_NumberBins + 1]);
    Kbari.push_back(X[i + 3*m_NumberBins + 1]);
  }
  if(m_FitDeltaKpi) {
    DeltaKpi = X[4*m_NumberBins + 1];
  }
  m_cisiLikelihood.PrintComparison(X[0], ci, si, Ki, Kbari, DeltaKpi);
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
  const std::vector<double> Generator_Ki = GetGeneratorcisi("K");
  const std::vector<double> Generator_Kbari = GetGeneratorcisi("Kbar");
  double Generator_DeltaKpi = m_Settings.getD("Generator_DeltaKpi");
  std::vector<double> ci_value(m_NumberBins), si_value(m_NumberBins);
  std::vector<double> ci_err(m_NumberBins), si_err(m_NumberBins);
  std::vector<double> ci_pull(m_NumberBins), si_pull(m_NumberBins);
  std::vector<double> Ki_value(m_NumberBins), Kbari_value(m_NumberBins);
  std::vector<double> Ki_err(m_NumberBins), Kbari_err(m_NumberBins);
  std::vector<double> Ki_pull(m_NumberBins), Kbari_pull(m_NumberBins);
  double DeltaKpi_value, DeltaKpi_err, DeltaKpi_pull;
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
    Tree.Branch(("K" + std::to_string(i + 1) + "_value").c_str(), &Ki_value[i]);
    Tree.Branch(("Kbar" + std::to_string(i + 1) + "_value").c_str(), &Kbari_value[i]);
    Tree.Branch(("K" + std::to_string(i + 1) + "_err").c_str(), &Ki_err[i]);
    Tree.Branch(("Kbar" + std::to_string(i + 1) + "_err").c_str(), &Kbari_err[i]);
    Tree.Branch(("K" + std::to_string(i + 1) + "_pull").c_str(), &Ki_pull[i]);
    Tree.Branch(("Kbar" + std::to_string(i + 1) + "_pull").c_str(), &Kbari_pull[i]);
  }
  if(m_FitDeltaKpi) {
    Tree.Branch("DeltaKpi_value", &DeltaKpi_value);
    Tree.Branch("DeltaKpi_err", &DeltaKpi_err);
    Tree.Branch("DeltaKpi_pull", &DeltaKpi_pull);
  }
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(-1);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const double BF_KKpipi = x[0];
    std::vector<double> ci, si, Ki, Kbari;
    double DeltaKpi = 0.0;
    for(std::size_t i = 0; i < m_NumberBins; i++) {
      ci.push_back(x[i + 0*m_NumberBins + 1]);
      si.push_back(x[i + 1*m_NumberBins + 1]);
      Ki.push_back(x[i + 2*m_NumberBins + 1]);
      Kbari.push_back(x[i + 3*m_NumberBins + 1]);
    }
    if(m_FitDeltaKpi) {
      DeltaKpi = x[4*m_NumberBins + 1];
    }
    return m_cisiLikelihood.CalculateToyLogLikelihood(BF_KKpipi,
						      ci, si,
						      Ki, Kbari,
						      DeltaKpi);
  };
  std::size_t StatsMultiplier = m_Settings.getI("StatsMultiplier");
  std::size_t NumberOfToysToSkip = m_Settings.getI("NumberOfToysToSkip");
  for(std::size_t i = 0; i < NumberOfToysToSkip; i++) {
    std::cout << "Skipped toy " << i << "\n";
    cisiLikelihoodRef.GenerateToy(Generator_BF_KKpipi,
				  Generator_ci,
				  Generator_si,
				  Generator_Ki,
				  Generator_Kbari,
				  Generator_DeltaKpi,
				  StatsMultiplier);
  }
  std::size_t NumberToys = m_Settings.getI("NumberToys");
  for(std::size_t i = 0; i < NumberToys; i++) {
    std::cout << "Toy " << i << "\n";
    cisiLikelihoodRef.GenerateToy(Generator_BF_KKpipi,
				  Generator_ci,
				  Generator_si,
				  Generator_Ki,
				  Generator_Kbari,
				  Generator_DeltaKpi,
				  StatsMultiplier);
    std::size_t NumberParameters = 4*m_NumberBins + 1;
    if(m_FitDeltaKpi) {
      NumberParameters++;
    }
    ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
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
      ci_value[i] = X[i + 0*m_NumberBins + 1];
      si_value[i] = X[i + 1*m_NumberBins + 1];
      Ki_value[i] = X[i + 2*m_NumberBins + 1];
      Kbari_value[i] = X[i + 3*m_NumberBins + 1];
      ci_err[i] = E[i + 0*m_NumberBins + 1];
      si_err[i] = E[i + 1*m_NumberBins + 1];
      Ki_err[i] = E[i + 2*m_NumberBins + 1];
      Kbari_err[i] = E[i + 3*m_NumberBins + 1];
      double ErrorLow, ErrorHigh;
      Minimiser.GetMinosError(i + 1 + 0*m_NumberBins, ErrorLow, ErrorHigh, 0);
      ci_pull[i] = ci_value[i] - Generator_ci[i];
      ci_pull[i] /= ci_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      //ci_pull[i] /= ci_err[i];
      Minimiser.GetMinosError(i + 1 + 1*m_NumberBins, ErrorLow, ErrorHigh, 0);
      si_pull[i] = si_value[i] - Generator_si[i];
      si_pull[i] /= si_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      //si_pull[i] /= si_err[i];
      Minimiser.GetMinosError(i + 1 + 2*m_NumberBins, ErrorLow, ErrorHigh, 0);
      Ki_pull[i] = Ki_value[i] - Generator_Ki[i];
      Ki_pull[i] /= Ki_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      //Ki_pull[i] /= Ki_err[i];
      Minimiser.GetMinosError(i + 1 + 3*m_NumberBins, ErrorLow, ErrorHigh, 0);
      Kbari_pull[i] = Kbari_value[i] - Generator_Kbari[i];
      Kbari_pull[i] /= Kbari_pull[i] <= 0.0 ? ErrorHigh : -ErrorLow;
      //Kbari_pull[i] /= Kbari_err[i];
    }
    if(m_FitDeltaKpi) {
      DeltaKpi_value = X[4*m_NumberBins + 1];
      DeltaKpi_err = E[4*m_NumberBins + 1];
      DeltaKpi_pull = (DeltaKpi_value - Generator_DeltaKpi)/DeltaKpi_err;
    }
    Tree.Fill();
  }
  Tree.Write();
  File.Close();
}

void cisiFitter::SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const {
  Minimiser.SetVariable(0, "BF_KKpipi", 0.00247, 0.1);
  Minimiser.SetVariableLimits(0, 0.001, 0.010);
  Minimiser.FixVariable(0);
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
  const double KMin = m_Settings.getD("K_min");
  const double KMax = m_Settings.getD("K_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "K" + std::to_string(i);
    Minimiser.SetVariable(i + 2*m_NumberBins, Name, 0.0, 1.0);
    Minimiser.SetVariableLimits(i + 2*m_NumberBins, KMin, KMax);
  }
  const double KbarMin = m_Settings.getD("Kbar_min");
  const double KbarMax = m_Settings.getD("Kbar_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "Kbar" + std::to_string(i);
    Minimiser.SetVariable(i + 3*m_NumberBins, Name, 0.0, 1.0);
    Minimiser.SetVariableLimits(i + 3*m_NumberBins, KbarMin, KbarMax);
  }
  if(m_FitDeltaKpi) {
    Minimiser.SetVariable(4*m_NumberBins + 1, "DeltaKpi", 0.0, 10.0);
    Minimiser.SetVariableLimits(4*m_NumberBins + 1, -180.0, 360.0);
  }
}

std::vector<double> cisiFitter::GetGeneratorcisi(const std::string &c_or_s_or_K) const {
  std::vector<double> cisi(m_NumberBins);
  const std::string Prefix = "Generator_" + c_or_s_or_K;
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
