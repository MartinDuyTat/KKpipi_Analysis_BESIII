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
#include"TGraphErrors.h"
#include"TGraph2D.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"Settings.h"
#include"cisiLikelihood.h"
#include"cisiFitter.h"
#include"Utilities.h"
#include"cisiFitterParameters.h"
#include"GammaLikelihood.h"
#include"GammaFitterParameters.h"

cisiFitter::cisiFitter(const Settings &settings):
  m_cisiLikelihood(settings),
  m_NumberBins(settings["BinningScheme"].getI("NumberBins")),
  m_Settings(settings),
  m_FitDeltaKpi(settings.getB("FitDeltaKpi")) {
}

void cisiFitter::Minimise() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(2);
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
  Minimiser.SetPrintLevel(2);
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

void cisiFitter::GenerateFeldmanCousinsYields(std::size_t Parameter) const {
  if(Parameter < 1 || Parameter > 4) {
    throw std::invalid_argument("Feldman Cousins scan parameters must be between 1 and 4");
  }
  if(!std::filesystem::exists("FeldmanCousinsYields") ||
     !std::filesystem::is_directory("FeldmanCousinsYields")) {
    std::filesystem::create_directory("FeldmanCousinsYields");
  }
  double Increment = m_Settings.getD("FeldmanCousins_Increment");
  int NumberIncrements = m_Settings.getI("FeldmanCousins_NumberIncrements");
  const auto GeneratorParameters = GetGeneratorValues();
  for(int i = -NumberIncrements; i <= NumberIncrements; i++) {
    if(i == 0) {
      continue;
    }
    auto NewGeneratorParameters = GeneratorParameters;
    NewGeneratorParameters.m_si[Parameter - 1] += Increment*i;
    std::string Filename = "FeldmanCousinsYields/ScanParameter";
    Filename += std::to_string(Parameter) + "_Increment";
    Filename += std::string(i > 0 ? "P" : "M") + std::to_string(TMath::Abs(i));
    Filename += "_GeneratorYields.txt";
    std::ofstream File(Filename);
    m_cisiLikelihood.SavePredictedBinYields(File, NewGeneratorParameters);
    File.close();
  }
}

void cisiFitter::FitFeldmanCousinsToy(const std::string &ToyName,
				      int ToyNumber,
				      std::size_t Parameter) const {
  int CovStatus = -1;
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(2);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const cisiFitterParameters Parameters(x, m_NumberBins);
    return m_cisiLikelihood.CalculateLogLikelihood(Parameters);
  };
  std::cout << "Feldman Cousins toy " << ToyNumber << "\n";
  // First do nominal fit to toy
  m_cisiLikelihood.LoadFeldmanCousinsDataset(ToyName, ToyNumber);
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
  // Save chi2 at minimum in the nominal fit
  const double MinValue = Minimiser.MinValue();
  // Stupid way of finding generator value of si
  std::string Name = ToyName;
  Name = Name.substr(Name.find('_') + 10, std::string::npos);
  Name = Name.substr(0, Name.find('_'));
  int Sign = Name[0] == 'P' ? +1 : -1;
  Name = Name.substr(1, std::string::npos);
  int Steps = std::stoi(Name);
  double Increment = m_Settings.getD("FeldmanCousins_Increment");
  double Gen_si = m_Settings.getD("Generator_s" + std::to_string(Parameter));
  double New_si = Gen_si + Sign*Steps*Increment;
  // Fix si value to scan point and repeat fit
  Minimiser.SetVariableValue(m_NumberBins + Parameter, New_si);
  Minimiser.FixVariable(m_NumberBins + Parameter);
  Counter = 0;
  CovStatus = -1;
  while(CovStatus != 3 && Counter < 5) {
    if(Counter != 0) {
      std::cout << "Fit did not converge, fitting again ";
      std::cout << "(" << Counter << ")\n";
    }
    Minimiser.Minimize();
    CovStatus = Minimiser.CovMatrixStatus();
    Counter++;
  }
  // Save new chi2 at the scan point and save the difference
  const double MinValue_FixedParam = Minimiser.MinValue();
  Minimiser.ReleaseVariable(m_NumberBins + Parameter);
  if(!std::filesystem::exists("FeldmanCousinsResults") ||
     !std::filesystem::is_directory("FeldmanCousinsResults")) {
    std::filesystem::create_directory("FeldmanCousinsResults");
  }
  const std::string Filename = "FeldmanCousinsResults/" + ToyName
                             + std::to_string(ToyNumber) + ".txt";
  std::ofstream OutputFile(Filename);
  OutputFile << MinValue_FixedParam - MinValue << "\n";
  OutputFile.close();
}

void cisiFitter::FeldmanCousinsDataScan() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(2);
  auto &cisiLikelihoodRef = m_cisiLikelihood;
  auto LikelihoodFunction = [=, &cisiLikelihoodRef] (const double *x) {
    const cisiFitterParameters Parameters(x, m_NumberBins);
    return m_cisiLikelihood.CalculateLogLikelihood(Parameters);
  };
  const std::size_t NumberParameters = 4*m_NumberBins + 3;
  ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  Minimiser.Minimize();
  const double GlobalMinimum = Minimiser.MinValue();
  std::size_t Parameter = m_Settings.getI("FeldmanCousins_Parameter");
  if(Parameter == 0 || Parameter > 4) {
    throw std::runtime_error("Parameter must be between 1 and 4");
  }
  const double Gen_si = m_Settings.getD("Generator_s" + std::to_string(Parameter));
  double Increment = m_Settings.getD("FeldmanCousins_Increment");
  int NumberIncrements = m_Settings.getI("FeldmanCousins_NumberIncrements");
  Minimiser.FixVariable(m_NumberBins + Parameter);
  std::string Filename = m_Settings.get("FeldmanCousinsDataScanFilename");
  Filename += "_" + std::to_string(Parameter) + ".txt";
  std::ofstream File(Filename);
  for(int i = -NumberIncrements; i <= NumberIncrements; i++) {
    if(i == 0) {
      continue;
    }
    const double New_si = Gen_si + i*Increment;
    Minimiser.SetVariableValue(m_NumberBins + Parameter, New_si);
    Minimiser.Minimize();
    const double NewMinimum = Minimiser.MinValue();
    File << New_si << " " << NewMinimum - GlobalMinimum << "\n";
  }
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
    std::vector<double> Ri, Rbari;
    for(int Bin = m_NumberBins; Bin >= 1; Bin--) {
      Values.push_back(m_Settings.getD("Generator_R_M" + std::to_string(Bin)));
    }
    for(int Bin = 1; Bin < static_cast<int>(m_NumberBins); Bin++) {
      Values.push_back(m_Settings.getD("Generator_R_P" + std::to_string(Bin)));
    }
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
  Results.GetYaxis()->SetTitleOffset(1.2);
  Results.SetTitle(";r_{D}^{K#pi}cos(#delta_{D}^{K#pi});r_{D}^{K#pi}sin(#delta_{D}^{K#pi})");
  Results.Draw("AP");
  // Draw HFLAV average
  double rDCosDeltaKpi_HFLAV = -0.05857;
  double rDCosDeltaKpi_err_HFLAV = 0.00017;
  double rDSinDeltaKpi_HFLAV = -0.0073;
  double rDSinDeltaKpi_err_HFLAV = 0.0093;
  TGraphErrors HFLAV(1, &rDCosDeltaKpi_HFLAV, &rDSinDeltaKpi_HFLAV,
		     &rDCosDeltaKpi_err_HFLAV, &rDSinDeltaKpi_err_HFLAV);
  HFLAV.SetMarkerStyle(20);
  HFLAV.SetLineWidth(3);
  HFLAV.Draw("PE SAME");
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

void cisiFitter::MinimiseWithGamma() const {
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(2);
  Minimiser.SetMaxFunctionCalls(100000);
  GammaLikelihood gammaLikelihood(m_Settings, m_cisiLikelihood);
  auto LikelihoodFunction = [=, &gammaLikelihood] (const double *x) {
    const GammaFitterParameters Parameters(x, m_NumberBins);
    return gammaLikelihood.CalculateLogLikelihood(Parameters);
  };
  // Number of parameters in the fit:
  // 6 CP observables for gamma
  // 2NBins - 1 R_i for LHCb
  // 4 yield normalisations
  // NBins ci
  // NBins si
  // 2*NBins - 1 R_i
  // 2 DeltaKpi
  // 1 BF of KKpipi
  // 1 BF of KKpipi for KLpipi tag
  const std::size_t NumberParameters = 6*m_NumberBins + 12;
  ROOT::Math::Functor fcn(LikelihoodFunction, NumberParameters);
  Minimiser.SetFunction(fcn);
  SetupMinimiser(Minimiser);
  /*Minimiser.SetVariable(4*m_NumberBins + 3, "xMinus", 0.08, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 3, -0.4, 0.4);
  Minimiser.SetVariable(4*m_NumberBins + 4, "yMinus", -0.03, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 4, -0.4, 0.4);
  Minimiser.SetVariable(4*m_NumberBins + 5, "xPlus", -0.13, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 5, -0.4, 0.4);
  Minimiser.SetVariable(4*m_NumberBins + 6, "yPlus", -0.04, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 6, -0.4, 0.4);
  Minimiser.SetVariable(4*m_NumberBins + 7, "xXi", -0.03, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 7, -0.4, 0.4);
  Minimiser.SetVariable(4*m_NumberBins + 8, "yXi", -0.02, 0.1);
  Minimiser.SetVariableLimits(4*m_NumberBins + 8, -0.4, 0.4);*/

  //Minimiser.SetVariable(4*m_NumberBins + 3, "gamma", 112.0, 10.0);
  Minimiser.SetVariable(4*m_NumberBins + 3, "gamma", 90.0, 20.0);
  Minimiser.SetVariableLimits(4*m_NumberBins + 3, 0.0, 180.0);
  //Minimiser.SetVariable(4*m_NumberBins + 4, "deltaB_dk", 84.0, 10.0);
  Minimiser.SetVariable(4*m_NumberBins + 4, "deltaB_dk", 90.0, 20.0);
  Minimiser.SetVariableLimits(4*m_NumberBins + 4, 0.0, 180.0);
  Minimiser.SetVariable(4*m_NumberBins + 5, "rB_dk", 0.1, 0.05);
  Minimiser.SetVariableLimits(4*m_NumberBins + 5, 0.0, 1.0);
  Minimiser.SetVariable(4*m_NumberBins + 6, "deltaB_dpi", 270.0, 50.0);
  //Minimiser.SetVariable(4*m_NumberBins + 6, "deltaB_dpi", 270.0, 10.0);
  Minimiser.SetVariableLimits(4*m_NumberBins + 6, 180.0, 360);
  Minimiser.SetVariable(4*m_NumberBins + 7, "rB_dpi", 0.005, 0.005);
  Minimiser.SetVariableLimits(4*m_NumberBins + 7, 0.0, 0.1);

  Minimiser.SetVariable(4*m_NumberBins + 8, "dummy", 0.0, 0.02);
  //Minimiser.SetVariableLimits(4*m_NumberBins + 8, -1.0, 1.0);
  Minimiser.FixVariable(4*m_NumberBins + 8);
  int Counter = -m_NumberBins;
  for(std::size_t i = 4*m_NumberBins + 9; i < 6*m_NumberBins + 8; i++) {
    std::string Name = Counter > 0 ? "P" : "M";
    Name += std::to_string(TMath::Abs(Counter++));
    if(Counter == 0) {
      Counter++;
    }
    Minimiser.SetVariable(i, ("R_LHCb_" + Name).c_str(), 0.5, 1.0);
    Minimiser.SetVariableLimits(i, 0.0, 1.0);
  }
  Minimiser.SetVariable(6*m_NumberBins + 8, "BMinusDKYield", 1500.0, 0.1);
  Minimiser.SetVariableLimits(6*m_NumberBins + 8, 100.0, 5000.0);
  Minimiser.SetVariable(6*m_NumberBins + 9, "BPlusDKYield", 1500.0, 0.1);
  Minimiser.SetVariableLimits(6*m_NumberBins + 9, 100.0, 5000.0);
  Minimiser.SetVariable(6*m_NumberBins + 10, "BMinusDpiYield", 20000.0, 0.1);
  Minimiser.SetVariableLimits(6*m_NumberBins + 10, 100.0, 70000.0);
  Minimiser.SetVariable(6*m_NumberBins + 11, "BPlusDpiYield", 20000.0, 0.1);
  Minimiser.SetVariableLimits(6*m_NumberBins + 11, 100.0, 70000.0);
  if(m_Settings.getB("cisiFixed")) {
    for(std::size_t i = 0; i < 4*m_NumberBins + 3; i++) {
      Minimiser.FixVariable(i);
    }
  }
  Minimiser.Minimize();
  /*for(std::size_t i = 4*m_NumberBins + 3; i <= 4*m_NumberBins + 8; i++) {
    for(std::size_t j = 4*m_NumberBins + 3; j <= 4*m_NumberBins + 8; j++) {
      std::cout << Minimiser.Correlation(i, j);
      if(j < 4*m_NumberBins + 8) {
	std::cout << ", ";
      } else {
	std::cout << ",\n";
      }
    }
  }*/
  // Create and draw contours
  TCanvas c("c", "", 1200, 900);
  std::vector<double> x, y, z;
  Minimiser.SetPrintLevel(-1);
  for(double gamma = 40.0; gamma <= 170.0; gamma += 10.0) {
    for(double deltaB_dk = 30.0; deltaB_dk <= 160.0; deltaB_dk += 10.0) {
      std::cout << "Scanning gamma = " << gamma << ", delta_B_dk = " << deltaB_dk << "\n";
      Minimiser.SetVariableValue(4*m_NumberBins + 3, gamma);
      Minimiser.SetVariableValue(4*m_NumberBins + 4, deltaB_dk);
      Minimiser.FixVariable(4*m_NumberBins + 3);
      Minimiser.FixVariable(4*m_NumberBins + 4);
      Minimiser.Minimize();
      x.push_back(gamma);
      y.push_back(deltaB_dk);
      z.push_back(Minimiser.MinValue());
    }
  }
  TGraph2D Graph2d(x.size(), x.data(), y.data(), z.data());
  Graph2d.SetTitle("Simultaneous LHCb and BESIII fit;#gamma;#delta_{B}^{DK}");
  Graph2d.Draw("COLZ");
  /*const std::size_t Points = 31;
  Minimiser.SetPrintLevel(-1);
  std::vector<TGraph> Contours;
  std::vector<double> ErrorDefs_Plot{2.30, 6.18, 11.83};
  for(double ErrorDef : ErrorDefs_Plot) {
    Minimiser.SetErrorDef(ErrorDef);
    for(std::size_t Index = 4*m_NumberBins + 3;
	Index <= 4*m_NumberBins + 4;
	Index += 2) {
      std::vector<double> x(Points), y(Points);
      unsigned int nPoints = Points - 1;
      Minimiser.Contour(Index, Index + 1, nPoints, x.data(), y.data());
      x.back() = x[0];
      y.back() = y[0];
      Contours.push_back(TGraph(Points, x.data(), y.data()));
      Contours.back().SetLineWidth(3);
    }
  }
  Contours[0].GetXaxis()->SetLimits(0.0, 180.0);
  Contours[0].GetYaxis()->SetRangeUser(0.0, 180.0);
  Contours[0].SetTitle("Simultaneous LHCb and BESIII fit;#gamma;#delta_{B}^{DK}");
  Contours[0].Draw("AL");
  for(auto &Contour : Contours) {
    Contour.Draw("L SAME");
  }*/
  c.SaveAs("Contours_gamma_deltaB.pdf");
  /*Minimiser.SetPrintLevel(-1);
  double ErrorPlus, ErrorMinus;
  std::vector<double> ErrorDefs_Minos{1.0, 4.0, 9.0};
  for(double ErrorDef : ErrorDefs_Minos) {
    Minimiser.SetErrorDef(ErrorDef);
    Minimiser.GetMinosError(4*m_NumberBins + 3, ErrorMinus, ErrorPlus);
    std::cout << "DeltaLL = " << ErrorDef << " gamma uncertainty: +";
    std::cout << ErrorPlus << ", " << ErrorMinus << "\n";
  }*/
}
