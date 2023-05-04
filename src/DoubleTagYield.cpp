// Martin Duy Tat 26th November 2021

#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include"TTree.h"
#include"TPad.h"
#include"TCanvas.h"
#include"TAxis.h"
#include"TLine.h"
#include"TCut.h"
#include"TFile.h"
#include"TRandom.h"
#include"TLatex.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooSimultaneous.h"
#include"RooPlot.h"
#include"RooHist.h"
#include"RooMsgService.h"
#include"RooHelpers.h"
#include"RooStats/SPlot.h"
#include"RooStats/RooStatsUtils.h"
#include"DoubleTagYield.h"
#include"Settings.h"
#include"BinnedDataLoader.h"
#include"BinnedFitModel.h"
#include"Category.h"
#include"Utilities.h"
#include"Bes3plotstyle.h"

using namespace RooFit;

DoubleTagYield::DoubleTagYield(const Settings &settings, TTree *Tree):
  m_SignalMBC(settings.get("FitVariable").c_str(), "",
	      settings.getD("FitRange_low"),
	      settings.getD("FitRange_high")),
  m_Settings(settings), m_Tree(Tree),
  m_DataLoader(m_Settings, m_Tree, &m_SignalMBC),
  m_FitModel(m_Settings, &m_SignalMBC) {
  for(int i = 0; i < 2; i++) {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Plotting);
  }
  m_SignalMBC.setBins(500, "cache");
}

void DoubleTagYield::DoFit() {
  RooDataSet *DataSet = m_DataLoader.GetDataSet();
  RooSimultaneous *Model = m_FitModel.GetPDF();
  std::vector<std::string> Categories =
    m_DataLoader.GetCategoryObject().GetCategories();
  // Perform an initial fit
  RooArgSet *Parameters = Model->getParameters(m_SignalMBC);
  m_InitialParameters = Parameters->snapshot();
  int nCPUs = 1;
  if(Categories.size() > 1) {
    nCPUs = 4;
  }
  auto Result = Model->fitTo(*DataSet, Save(), NumCPU(nCPUs),
			     Strategy(2), Minos(m_FitModel.m_SignalYields),
			     Minimizer("Minuit2","migrad"));
  // Perform a second fit if necessary
  if(m_Settings.contains("SecondFit") && m_Settings.getB("SecondFit")) {
    Result = Model->fitTo(*DataSet, Save(), NumCPU(nCPUs),
			  Strategy(2), Minos(m_FitModel.m_SignalYields),
			  Minimizer("Minuit2","migrad"));
  }
  Result->Print("V");
  m_ParametersAfterFit = Parameters->snapshot();
  PlotProjections();
  SaveSignalYields(m_Settings.get("FittedSignalYieldsFile"), Result);
  // Smear peaking backgrounds for systematics studies
  if(m_Settings.getB("YieldSystematics")) {
    DoSystematicsFits();
  }
  // Perform toy fits with generator signal yields
  if(m_Settings.getB("ToyFits")) {
    DoToyFits();
  }
  // sPlot
  if(m_Settings.contains("sPlotReweight") &&
     m_Settings.getB("sPlotReweight")) {
    sPlotReweight(*DataSet);
  }
}

void DoubleTagYield::DoToyFits() {
  RooSimultaneous *Model = m_FitModel.GetPDF();
  RooArgSet *Parameters = Model->getParameters(m_SignalMBC);
  auto CatObject = m_DataLoader.GetCategoryObject();
  int nCPUs = 1;
  if(CatObject.GetCategories().size() > 1) {
    nCPUs = 4;
  }
  int NumberToys = m_Settings.getI("NumberToys");
  RooHelpers::LocalChangeMsgLevel changeMsgLvl(RooFit::ERROR);
  for(int i = 0; i < NumberToys; i++) {
    std::cout << "Running toy " << i << "...\n";
    *Parameters = *m_ParametersAfterFit;
    m_FitModel.SetGeneratorYields();
    gRandom->SetSeed(m_Settings.getI("Seed"));
    auto ToyDataset = Model->generate(RooArgSet(m_SignalMBC,
						*CatObject.GetCategoryVariable()),
				      Extended());
    auto Result = Model->fitTo(*ToyDataset, Save(), NumCPU(nCPUs),
			       Strategy(2), Minos(true),
			       PrintLevel(-1),
			       Minimizer("Minuit2","migrad"));
    std::string Filename = "ToyFitResults/Toy";
    Filename += std::to_string(i) + ".txt";
    SaveSignalYields(Filename, Result);
    std::cout << "Toy " << i << " done!\n";
  }
}

void DoubleTagYield::DoSystematicsFits() {
  std::ofstream OutputFile(m_Settings.get("FittedSignalYieldsFile"),
			   std::ios_base::app);
  OutputFile << "\n* Systematic uncertainties\n\n";
  std::vector<std::string> Categories =
    m_DataLoader.GetCategoryObject().GetCategories();
  std::map<std::string, double> SystError;
  TMatrixT<double> SystCovMatrix(Categories.size(), Categories.size());
  int PeakingBackgrounds =
    m_Settings["MBC_Shape"].getI(m_Settings.get("Mode") + "_PeakingBackgrounds");
  if(PeakingBackgrounds > 0) {
    std::map<std::string, std::vector<double>> FittedYields;
    for(const auto &Category : Categories) {
      FittedYields.insert({Category, std::vector<double>()});
    }
    gRandom->SetSeed(m_Settings.getI("Seed"));
    m_FitModel.PrepareSmearing();
    int SuccessfulFits = 0;
    RooDataSet *DataSet = m_DataLoader.GetDataSet();
    RooSimultaneous *Model = m_FitModel.GetPDF();
    RooArgSet *Parameters = Model->getParameters(m_SignalMBC);
    int nCPUs = 1;
    if(Categories.size() > 1) {
      nCPUs = 4;
    }
    for(int i = 0; i < m_Settings.getI("NumberRuns"); i++) {
      std::cout << "Starting systematics fit number: " << i << "\n";
      *Parameters = *m_InitialParameters;
      m_FitModel.SmearPeakingBackgrounds();
      auto Result = Model->fitTo(*DataSet, Strategy(2), Save(), NumCPU(nCPUs));
      Result->Print("V");
      if(Result->status() == 0 && Result->covQual() == 3) {
	SuccessfulFits++;
      } else {
	continue;
      }
      for(const auto &Category : Categories) {
	double Yield = m_FitModel.m_Yields[Category + "_SignalYield"]->getVal();
	FittedYields[Category].push_back(Yield);
      }
    }
    std::cout << "Number of successfull fits: " << SuccessfulFits << "\n";
    for(const auto &Category : Categories) {
      SystError[Category] = TMath::RMS(FittedYields[Category].begin(),
				       FittedYields[Category].end());
    }
    auto CategoryObj = m_DataLoader.GetCategoryObject();
    for(const auto &Category_x : Categories) {
      for(const auto &Category_y : Categories) {
	double Covariance = Utilities::Covariance(FittedYields[Category_x],
						  FittedYields[Category_y]);
	int index_x = CategoryObj.GetCategoryIndex(Category_x);
	int index_y = CategoryObj.GetCategoryIndex(Category_y);
	SystCovMatrix(index_x, index_y) = Covariance;
      }
    }
  } else {
    for(const auto &Category : Categories) {
      SystError[Category] = 0.0;
    }
  }
  for(const auto &Category : Categories) {
    std::string YieldName = Category + "_SignalYield_PeakingBackgrounds";
    OutputFile << YieldName << "_syst_err " << SystError[Category] << "\n";
  }
  OutputFile.close();
  if(Categories.size() > 1) {
    TFile SystCovMatrixFile("PeakingBackground_CovMatrix.root", "RECREATE");
    SystCovMatrixFile.cd();
    SystCovMatrix.Write("CovMatrix");
    SystCovMatrixFile.Close();
  }
}

void DoubleTagYield::PlotProjections() {
  using namespace RooFit;
  SetStyle();
  SetPrelimStyle();
  auto Model = m_FitModel.GetPDF();
  auto category = m_DataLoader.GetCategoryObject();
  auto CategoryVariable = category.GetCategoryVariable();
  RooDataSet *DataSet = m_DataLoader.GetDataSet();
  for(const auto &Category : category.GetCategories()) {
    int SignalBin = category.GetSignalBinNumber(Category);
    int TagBin = category.GetTagBinNumber(Category);
    TCanvas c1((Category + "_c1").c_str(), "", 1600, 1600);
    TPad Pad1((Category + "_Pad1").c_str(), "", 0.0, 0.25, 1.0, 1.0);
    TPad Pad2((Category + "_Pad2").c_str(), "", 0.0, 0.0, 1.0, 0.25);
    Pad1.Draw();
    Pad2.Draw();
    Pad1.cd();
    RooPlot *Frame = m_SignalMBC.frame();
    FormatAxis(Frame->GetXaxis());
    FormatAxis(Frame->GetYaxis());
    if(m_Settings.contains("No_x_axis_tick_label") &&
       m_Settings.getB("No_x_axis_tick_label")) {
      Frame->GetXaxis()->SetLabelSize(0);
      std::cout << "Removing tick labels\n";
    }
    std::string TagMode = m_Settings.get("Mode");
    TLatex Text;
    Text.SetTextFont(42);
    Text.SetTextSize(0.09);
    Text.SetTextColor(kBlack);
    Text.SetNDC(true);
    std::string LabelText = Utilities::GetTagNameLaTeX(TagMode);
    std::string Title;
    if(SignalBin != 0) {
      LabelText += ", bin " + std::to_string(SignalBin);
    }
    if(TagMode == "KSpipi" || TagMode == "KSpipiPartReco" ||
       TagMode == "KLpipi" || TagMode == "KSKK" ||
       TagMode == "KLKK" || TagMode == "KKpipi") {
      LabelText = "#splitline{" + LabelText;
      LabelText += "}{Bin " + std::to_string(TagBin) + "}";
    }
    if(TagMode.substr(0, 2) == "KL" ||
       TagMode.find("PartReco") != std::string::npos) {
      Title += ";M_{ miss}^{ 2} (GeV^{2}/#it{c}^{4}); Events / ";
    } else if(TagMode == "KeNu") {
      Title += ";U_{ miss} (GeV/#it{c}^{2}); Events / ";
    } else {
      Title += ";M_{ BC} (GeV/#it{c}^{2}); Events / ";
    }
    Text.SetText(0.2, 0.8, LabelText.c_str());
    if(TagMode.substr(0, 2) == "KL" ||
       TagMode.find("PartReco") != std::string::npos) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(3);
      double RangePerBin = m_SignalMBC.getMax() - m_SignalMBC.getMin();
      RangePerBin /= m_Settings.getI("Bins_in_plots");
      ss << RangePerBin;
      Title += ss.str();
      Title += " GeV^{2}/#it{c}^{4}";
    } else {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(1);
      double RangePerBin = m_SignalMBC.getMax() - m_SignalMBC.getMin();
      RangePerBin *= 1000.0/m_Settings.getI("Bins_in_plots");
      ss << RangePerBin;
      Title += ss.str();
      Title += " MeV/#it{c}^{2}";
    }
    Frame->SetTitle(Title.c_str());
    std::string CategoryCut = std::string(CategoryVariable->GetName());
    CategoryCut += "==" + std::string(CategoryVariable->GetName());
    CategoryCut += "::" + Category;
    RooPlot *Data_RooPlot =
      DataSet->plotOn(Frame,
		      Binning(m_Settings.getI("Bins_in_plots")),
		      MarkerSize(3),
		      LineWidth(3),
		      Cut(CategoryCut.c_str()));
    auto Data_RooHist = Data_RooPlot->getHist();
    // Against my wishes, I had to remove data points of empty bins
    if(m_Settings.contains("HideEmptyBins") && m_Settings.getB("HideEmptyBins")) {
      for(int i = 0; i < Data_RooHist->GetN(); i++) {
	if(Data_RooHist->GetPointY(i) == 0.0) {
	  Data_RooHist->SetPointY(i, -1000.0);
	}
      }
    }
    FormatData(Data_RooHist);
    Data_RooHist->SetMinimum(0.0);
    if(TagMode == "KSpipiPartReco") {
      Frame->SetNdivisions(-405);
    }
    Model->plotOn(Frame,
		  LineColor(kRed),
		  LineWidth(3),
		  Slice(*CategoryVariable, Category.c_str()),
		  ProjWData(*CategoryVariable, *DataSet));
    RooHist *Pull = Frame->pullHist();
    if(m_FitModel.m_PeakingBackgroundShapes.size() > 0) {
      std::string PeakingList = "Combinatorial*,";
      for(auto iter = m_FitModel.m_PeakingBackgroundShapes.begin();
	  iter != m_FitModel.m_PeakingBackgroundShapes.end();
	  iter++) {
	if(iter != m_FitModel.m_PeakingBackgroundShapes.begin()) {
	  PeakingList += std::string(",");
	}
	PeakingList += iter->second->GetPDF()->GetName();
      }
      Model->plotOn(Frame,
		    FillStyle(1001),
		    LineColor(kGreen + 2),
		    FillColor(kGreen + 2),
		    LineWidth(3),
		    DrawOption("F"),
		    Slice(*CategoryVariable, Category.c_str()),
		    ProjWData(*CategoryVariable, *DataSet),
		    Components(PeakingList.c_str()));
    }
    Model->plotOn(Frame,
		  FillStyle(1001),
		  LineColor(kAzure + 6),
		  FillColor(kAzure + 6),
		  LineWidth(3),
		  DrawOption("F"),
		  Components("Combinatorial*"),
		  LineStyle(kDashed),
		  Slice(*CategoryVariable, Category.c_str()),
		  ProjWData(*CategoryVariable, *DataSet));
    Data_RooPlot =
      DataSet->plotOn(Frame,
		      Binning(m_Settings.getI("Bins_in_plots")),
		      MarkerSize(3),
		      LineWidth(3),
		      Cut(CategoryCut.c_str()));
    Data_RooHist = Data_RooPlot->getHist();
    // Against my wishes, I had to remove data points of empty bins
    if(m_Settings.contains("HideEmptyBins") && m_Settings.getB("HideEmptyBins")) {
      for(int i = 0; i < Data_RooHist->GetN(); i++) {
	if(Data_RooHist->GetPointY(i) == 0.0) {
	  Data_RooHist->SetPointY(i, -1000.0);
	}
      }
    }
    Frame->Draw();
    Text.Draw("SAME");
    Pad2.cd();
    RooPlot *PullFrame = m_SignalMBC.frame();
    PullFrame->addObject(Pull, "P");
    TLine *Line = new TLine(m_SignalMBC.getMin(), 0.0, m_SignalMBC.getMax(), 0.0);
    PullFrame->addObject(Line);
    PullFrame->SetMinimum(-5);
    PullFrame->SetMaximum(5);
    PullFrame->SetTitle(";;");
    PullFrame->Draw();
    std::string PlotFilename = m_Settings.get("MBCPlotFilenamePrefix");
    PlotFilename += "_" + Category + ".pdf";
    Pad1.SetFrameLineWidth(3);
    c1.Draw();
    c1.SaveAs(PlotFilename.c_str());
  }
}

void DoubleTagYield::SaveSignalYields(const std::string &Filename,
				      RooFitResult *Result) const {
  // Open output file
  std::ofstream Outfile(Filename);
  Outfile << std::fixed << std::setprecision(4);
  // Save fit status
  Outfile << "* KKpipi vs " << m_Settings.get("Mode");
  Outfile << " double tag yield fit results\n\n";
  Outfile << "status " << Result->status() << "\n";
  Outfile << "covQual " << Result->covQual() << "\n\n";
  if(m_Settings.getB("FullyReconstructed")) {
    Outfile << "* These yields are after the sideband has been";
    Outfile << " subtracted off the signal yield\n\n";
  }
  // Loop over all categories and save signal yields
  auto category = m_DataLoader.GetCategoryObject();
  for(const auto & cat : category.GetCategories()) {
    std::string Name = cat + "_SignalYield";
    double Sideband = 0.0;
    if(m_Settings.getB("FullyReconstructed")) {
      Sideband += GetSidebandYield(category.GetSignalBinNumber(cat),
				   category.GetTagBinNumber(cat));
    }
    auto YieldVariable = static_cast<RooRealVar*>(m_FitModel.m_Yields.at(Name));
    Outfile << Name << "          " << std::setw(8) << std::right;
    Outfile << YieldVariable->getVal() - Sideband << "\n";
    Outfile << Name << "_err      " << std::setw(8) << std::right;
    Outfile << YieldVariable->getError() << "\n";
    Outfile << Name << "_low_err  " << std::setw(8) << std::right;
    Outfile << YieldVariable->getErrorLo() << "\n";
    Outfile << Name << "_high_err " << std::setw(8) << std::right;
    Outfile << YieldVariable->getErrorHi() << "\n";
    if(m_Settings.getB("FullyReconstructed")) {
      Outfile << Name << "_sideband " << std::setw(8) << std::right;
      Outfile << Sideband << "\n";
    }
  }
  // Save other parameters as well
  Outfile << std::fixed << std::setprecision(6);
  for(const auto &FitParameter : m_FitModel.m_Parameters) {
    auto Parameter = FitParameter.second;
    const std::string Name = Parameter->GetName();
    Outfile << Name << "          " << std::setw(12) << std::right;
    Outfile << Parameter->getVal() << "\n";
    Outfile << Name << "_err      " << std::setw(12) << std::right;
    Outfile << Parameter->getError() << "\n";
    Outfile << Name << "_low_err  " << std::setw(12) << std::right;
    Outfile << Parameter->getErrorLo() << "\n";
    Outfile << Name << "_high_err " << std::setw(12) << std::right;
    Outfile << Parameter->getErrorHi() << "\n";
  }
  Outfile.close();
  // Loop over all categories and save correlation and covarience matrices
  std::size_t Size = category.GetCategories().size();
  if(Size > 0) {
    TFile File("RawYieldsCorrelationMatrix.root", "RECREATE");
    TMatrixTSym<double> CorrelationMatrix(Size);
    TMatrixTSym<double> CovarianceMatrix(Size);
    const auto categories = category.GetCategories();
    for(std::size_t i = 0; i < Size; i++) {
      const std::string Name_i = categories[i] + "_SignalYield";
      auto YieldVariable_i =
	static_cast<RooRealVar*>(m_FitModel.m_Yields.at(Name_i));
      const double Error_i = YieldVariable_i->getError();
      for(std::size_t j = 0; j < Size; j++) {
	const std::string Name_j = categories[j] + "_SignalYield";
	auto YieldVariable_j =
	  static_cast<RooRealVar*>(m_FitModel.m_Yields.at(Name_j));
	const double Error_j = YieldVariable_j->getError();
	CorrelationMatrix(i, j) = Result->correlation(Name_i.c_str(),
						      Name_j.c_str());
	CovarianceMatrix(i, j) = CorrelationMatrix(i, j)*Error_i*Error_j;
      }
    }
    CorrelationMatrix.Write("CorrelationMatrix");
    CovarianceMatrix.Write("CovarianceMatrix");
    File.Close();
  }
}

double DoubleTagYield::GetSidebandYield(int SignalBin, int TagBin) const {
  TCut SidebandCut("TagMBC > 1.84 && TagMBC < 1.85");
  SidebandCut = SidebandCut && TCut("SignalMBC > 1.86 && SignalMBC < 1.87");
  if(SignalBin != 0) {
    SidebandCut = SidebandCut &&
                  TCut(("SignalBin == " + std::to_string(SignalBin)).c_str());
  }
  if(TagBin != 0) {
    SidebandCut = SidebandCut &&
                  TCut(("TagBin == " + std::to_string(TagBin)).c_str());
  }
  return static_cast<double>(m_Tree->GetEntries(SidebandCut));
}

void DoubleTagYield::sPlotReweight(RooDataSet &Data) {
  auto FloatingParameters = m_FitModel.GetPDF()->getParameters(Data);
  for(auto Parameter : *FloatingParameters) {
    auto CastedParameter = dynamic_cast<RooRealVar*>(Parameter);
    if(CastedParameter) {
      CastedParameter->setConstant(true);
    }
  }
  RooArgList YieldParameters;
  for(const auto &CategoryString : m_FitModel.m_Category.GetCategories()) {
    std::string YieldName = CategoryString + "_SignalYield";
    auto YieldParameter =
      static_cast<RooRealVar*>(m_FitModel.m_Yields.at(YieldName));
    YieldParameter->setConstant(false);
    YieldParameters.add(*YieldParameter);
  }
  RooStats::SPlot("sData", "", Data, m_FitModel.GetPDF(), YieldParameters);
  TFile Outfile(m_Settings.get("sPlotFilename").c_str(), "RECREATE");
  auto Tree = RooStats::GetAsTTree(m_Settings.get("TreeName").c_str(),
				   m_Settings.get("TreeName").c_str(),
				   Data);
  Tree->Write();
  Outfile.Close();
}
