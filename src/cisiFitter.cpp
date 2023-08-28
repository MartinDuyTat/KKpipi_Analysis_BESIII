// Martin Duy Tat 13th October 2022

#include<fstream>
#include<filesystem>
#include<vector>
#include"TFile.h"
#include"TTree.h"
#include"TRandom.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"
#include"TLatex.h"
#include"TLine.h"
#include"TParameter.h"
#include"TEllipse.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"Settings.h"
#include"cisiLikelihood.h"
#include"cisiFitter.h"
#include"Utilities.h"
#include"cisiFitterParameters.h"

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
    const cisiFitterParameters Parameters(x, m_NumberBins);
    return m_cisiLikelihood.CalculateLogLikelihood(Parameters);
  };
  // Number of parameters in the fit:
  // NBins ci
  // NBins si
  // 2*NBins - 1 R_i
  // 2 DeltaKpi
  // 1 BF of KKpipi
  // 1 BF of KKpipi for KLpipi tag
  const std::size_t NumberParameters = 4*m_NumberBins + 3;
  ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  Minimiser.Minimize();
  const cisiFitterParameters Parameters(Minimiser.X(), m_NumberBins);
  m_cisiLikelihood.PrintComparison(Parameters);
  if(m_Settings.getB("PlotContours")) {
    Plot_cisi(Minimiser, Parameters.m_ci, Parameters.m_si);
  }
  if(m_FitDeltaKpi && m_Settings.getB("PlotDeltaKpiContour")) {
    Plot_DeltaKpi(Minimiser,
		  Parameters.m_rDcosDeltaKpi,
		  Parameters.m_rDsinDeltaKpi);
  }
  std::string ResultsFile;
  if(m_Settings.get("Systematics") == "None") {
    ResultsFile = m_Settings.get("FitResultsFile");
  } else {
    std::string SystDir = m_Settings.get("Systematics") + "Systematics";
    if(!std::filesystem::exists(SystDir) ||
       !std::filesystem::is_directory(SystDir)) {
      std::filesystem::create_directory(SystDir);
    }
    int FitNumber;
    if(m_Settings.get("Systematics") == "PeakingBackgrounds") {
      FitNumber = m_Settings.getI("FitNumber");
    } else {
      FitNumber = m_Settings.getI("SystSeed");
    }
    ResultsFile = SystDir + "/Fit" + std::to_string(FitNumber) + ".root";
  }
  SaveFitResults(Minimiser, ResultsFile);
}

void cisiFitter::RunToy(int ToyNumber) const {
  int CovStatus = -1;
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(4);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const cisiFitterParameters Parameters(x, m_NumberBins);
    return m_cisiLikelihood.CalculateLogLikelihood(Parameters);
  };
  std::cout << "Toy " << ToyNumber << "\n";
  m_cisiLikelihood.LoadToyDataset(static_cast<int>(ToyNumber));
  const std::size_t NumberParameters = 4*m_NumberBins + 3;
  ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  std::size_t Counter = 0;
  while(CovStatus != 3 && Counter < 5) {
    if(Counter != 0) {
      std::cout << "Fit did not converge, fitting again ";
      std::cout << "(" << Counter << ")\n";
    }
    Minimiser.Minimize();
    CovStatus = Minimiser.CovMatrixStatus();
    Counter++;
  }
  if(!std::filesystem::exists("ToyFitResults") ||
     !std::filesystem::is_directory("ToyFitResults")) {
    std::filesystem::create_directory("ToyFitResults");
  }
  const std::string Filename = "ToyFitResults/Toy"
                             + std::to_string(ToyNumber) + ".root";
  SaveFitResults(Minimiser, Filename);
}

void cisiFitter::SavePredictedYields() const {
  const std::string Filename = m_Settings.get("PredictedYieldsFile");
  std::ofstream File(Filename);
  const auto GeneratorParameters = GetGeneratorValues();
  m_cisiLikelihood.SavePredictedBinYields(File, GeneratorParameters);
  File.close();
}

void cisiFitter::SetupMinimiser(ROOT::Minuit2::Minuit2Minimizer &Minimiser) const {
  Minimiser.SetVariable(0, "BF_KKpipi", 0.00247, 0.1);
  Minimiser.SetVariableLimits(0, 0.001, 0.010);
  const double cMin = m_Settings.getD("c_min");
  const double cMax = m_Settings.getD("c_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "c" + std::to_string(i);
    const double Initial = 0.0;
    Minimiser.SetVariable(i, Name, Initial, 1.0);
    Minimiser.SetVariableLimits(i, cMin, cMax);
  }
  const double sMin = m_Settings.getD("s_min");
  const double sMax = m_Settings.getD("s_max");
  for(std::size_t i = 1; i <= m_NumberBins; i++) {
    const std::string Name = "s" + std::to_string(i);
    const double Initial = 0.0;
    Minimiser.SetVariable(i + m_NumberBins, Name, Initial, 1.0);
    Minimiser.SetVariableLimits(i + m_NumberBins, sMin, sMax);
    if(m_Settings.getB("Fix_si")) {
      Minimiser.FixVariable(i + m_NumberBins);
    }
  }
  std::size_t Counter = 2*m_NumberBins + 1;
  for(std::size_t i = m_NumberBins; i > 0; i--) {
    const std::string Name = "R_M" + std::to_string(i);
    Minimiser.SetVariable(Counter, Name, 0.0, 1.0);
    Minimiser.SetVariableLimits(Counter, 0.0, 1.0);
    Counter++;
  }
  for(std::size_t i = 1; i < m_NumberBins; i++) {
    const std::string Name = "R_P" + std::to_string(i);
    Minimiser.SetVariable(Counter, Name, 0.0, 1.0);
    Minimiser.SetVariableLimits(Counter, 0.0, 1.0);
    Counter++;
  }
  Minimiser.SetVariable(4*m_NumberBins, "BF_KKpipi_KLpipi", 0.00247, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins, 0.001, 0.010);
  if(m_Settings.get("TagModes").find("KLpipi") == std::string::npos) {
    Minimiser.FixVariable(4*m_NumberBins);
  }
  Minimiser.SetVariable(4*m_NumberBins + 1, "rDcosDeltaKpi", 0.0, 0.05);
  Minimiser.SetVariableLimits(4*m_NumberBins + 1, -0.3, 0.3);
  Minimiser.SetVariable(4*m_NumberBins + 2, "rDsinDeltaKpi", 0.0, 0.05);
  Minimiser.SetVariableLimits(4*m_NumberBins + 2, -0.5, 0.5);
  if(!m_FitDeltaKpi) {
    Minimiser.FixVariable(4*m_NumberBins + 1);
    Minimiser.FixVariable(4*m_NumberBins + 2);
  }
}

cisiFitterParameters cisiFitter::GetGeneratorValues() const {
  std::vector<double> Values;
  Values.push_back(m_Settings.getD("Generator_BF_KKpipi"));
  for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
    Values.push_back(m_Settings.getD("Generator_c" + std::to_string(Bin)));
  }
  for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
    Values.push_back(m_Settings.getD("Generator_s" + std::to_string(Bin)));
  }
  {
    std::vector<double> Ki, Kbari;
    for(std::size_t Bin = 1; Bin <= m_NumberBins; Bin++) {
      Ki.push_back(m_Settings.getD("Generator_K" + std::to_string(Bin)));
      Kbari.push_back(m_Settings.getD("Generator_Kbar" + std::to_string(Bin)));
    }
    const auto Ri = Utilities::ConvertKiToRi(Ki, Kbari);
    Values.insert(Values.end(), Ri.begin(), Ri.end() - 1);
  }
  Values.push_back(m_Settings.getD("Generator_BF_KKpipi"));
  {
    const double DeltaKpi = m_Settings.getD("Generator_DeltaKpi");
    const double rDKpi = m_Settings.getD("Generator_rDKpi");
    Values.push_back(rDKpi*TMath::Cos(DeltaKpi*TMath::Pi()/180.0));
    Values.push_back(rDKpi*TMath::Sin(DeltaKpi*TMath::Pi()/180.0));
  }
  return cisiFitterParameters(Values.data(), m_NumberBins);
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
  std::vector<double> ErrorDefs{2.30};
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

void cisiFitter::Plot_DeltaKpi(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
			       double rDcosDeltaKpi,
			       double rDsinDeltaKpi) const {
  // Make canvas
  TCanvas c("DeltaKpi_c", "", 900, 900);
  // Draw fit results
  TGraph Results(1, &rDcosDeltaKpi, &rDsinDeltaKpi);
  const double Boundary = m_Settings.getD("DeltaKpiPlotBoundary");
  Results.GetXaxis()->SetLimits(-Boundary, Boundary);
  Results.GetYaxis()->SetRangeUser(-Boundary, Boundary);
  Results.SetMarkerStyle(8);
  Results.SetTitle(";r_{D}cos(#delta_{K#pi});r_{D}sin(#delta_{K#pi})");
  Results.Draw("AP");
  // Draw circle representing rD
  TEllipse rDContour(0.0, 0.0, 0.05867, 0.05867);
  //rDContour.SetFillColor(17);
  rDContour.SetLineStyle(kDashed);
  rDContour.SetLineWidth(3);
  rDContour.Draw("SAME");
  // Create and draw contour
  const std::size_t Points = 101;
  Minimiser.SetPrintLevel(-1);
  std::vector<double> ErrorDefs{2.30, 6.18};
  TGraph *Contour = nullptr;
  for(double ErrorDef : ErrorDefs) {
    Minimiser.SetErrorDef(ErrorDef);
    std::vector<double> x(Points), y(Points);
    unsigned int nPoints = Points - 1;
    Minimiser.Contour(4*m_NumberBins + 1, 4*m_NumberBins + 2, nPoints,
		      x.data(), y.data());
    x.back() = x[0];
    y.back() = y[0];
    Contour = new TGraph(Points, x.data(), y.data());
    Contour->SetLineWidth(3);
    Contour->Draw("L SAME");
  }
  c.SaveAs("Contour_DeltaKpi.pdf");
  delete Contour;
}

void cisiFitter::SaveFitResults(ROOT::Minuit2::Minuit2Minimizer &Minimiser,
				const std::string &Filename) const {
  std::cout << "Saving fit results...\n";
  std::size_t NParameters = Minimiser.NDim();
  TFile File(Filename.c_str(), "RECREATE");
  {
    TParameter<int> Status("Status", Minimiser.Status());
    File.WriteObject(&Status, "Status");
    TParameter<int> CovStatus("CovStatus", Minimiser.CovMatrixStatus());
    File.WriteObject(&CovStatus, "CovStatus");
    TParameter<double> LL("LL", Minimiser.MinValue());
    File.WriteObject(&LL, "LL");
    std::cout << "LL and fit status saved\n";
  }
  {
    std::vector<std::string> VariableNames;
    for(std::size_t i = 0; i < NParameters; i++) {
      VariableNames.push_back(Minimiser.VariableName(i));
    }
    File.WriteObject(&VariableNames, "VariableNames");
    std::cout << "Variable names saved\n";
  }
  {
    const double *Values = Minimiser.X();
    const std::vector<double> FitValues(Values, Values + NParameters);
    File.WriteObject(&FitValues, "FitValues");
    std::cout << "Fit values saved\n";
  }
  {
    const double *Errors = Minimiser.Errors();
    const std::vector<double> FitUncertainties(Errors, Errors + NParameters);
    File.WriteObject(&FitUncertainties, "FitUncertainties");
    std::cout << "Fit uncertainties saved\n";
  }
  if(m_Settings.getB("RunMinos")) {
    Minimiser.SetPrintLevel(-1);
    std::vector<double> PlusUncertainties(NParameters),
                        MinusUncertainties(NParameters);
    for(std::size_t i = 0; i < NParameters; i++) {
      Minimiser.GetMinosError(i, MinusUncertainties[i], PlusUncertainties[i]);
    }
    File.WriteObject(&PlusUncertainties, "PlusUncertainties");
    File.WriteObject(&MinusUncertainties, "MinusUncertainties");
    std::cout << "Minos uncertainties saved\n";
  }
  {
    std::vector<double> CovarianceMatrix, CorrelationMatrix;
    for(std::size_t i = 0; i < NParameters; i++) {
      for(std::size_t j = 0; j < NParameters; j++) {
	CovarianceMatrix.push_back(Minimiser.CovMatrix(i, j));
	CorrelationMatrix.push_back(Minimiser.Correlation(i, j));
      }
    }
    File.WriteObject(&CovarianceMatrix, "CovMatrix");
    File.WriteObject(&CorrelationMatrix, "CorrMatrix");
    std::cout << "Covariance and correlation matrices saved\n";
  }
  File.Close();
  std::cout << "Fit results saved to " << Filename << "\n";
}
