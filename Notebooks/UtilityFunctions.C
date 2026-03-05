/**
 * Get branching fraction for a D decay
 * Mode Decay mode
 */
std::pair<double, double> GetBranchingFraction(const std::string &Mode) {
  std::string Filename = "/data/bes3/tat/KKpipi_StrongPhase_Analysis_4Bins/";
  Filename += "CommonInputs/PDG/BranchingFractions.txt";
  std::ifstream BFFile(Filename);
  std::string Line;
  std::pair<double, double> BF;
  while(std::getline(BFFile, Line)) {
    if(Line.empty()) {
      continue;
    }
    std::string Word;
    double Value;
    std::stringstream ss(Line);
    ss >> Word >> Value;
    if(Word == Mode) {
      BF.first = Value;
    } else if(Word == Mode + "_err") {
      BF.second = Value;
    }
  }
  return BF;
}

/**
 * Get F+ for a D decay
 * Mode Decay mode
 */
std::pair<double, double> GetFPlus(const std::string &Mode) {
  std::string Filename = "/data/bes3/tat/KKpipi_StrongPhase_Analysis_4Bins/";
  Filename += "CommonInputs/FPlus/FPlus_TagModes.txt";
  std::ifstream BFFile(Filename);
  std::string Line;
  std::pair<double, double> FPlus(-1.0, 0.0);
  while(std::getline(BFFile, Line)) {
    if(Line.empty()) {
      continue;
    }
    std::string Word;
    double Value;
    std::stringstream ss(Line);
    ss >> Word >> Value;
    if(Word == Mode) {
      FPlus.first = Value;
    } else if(Word == Mode + "_err") {
      FPlus.second = Value;
    }
  }
  return FPlus;
}

/**
 * Sum the event weights of a TTree
 * @param Tree The TTree
 * @param WeightName1 The name of first set of weights
 * @param WeightName2 The name of second set of weights
 * @param Weight1Fraction Fraction of weight 1
 * @param Cut Cut to apply, or leave empty
 */
double SumWeights(TTree *Tree,
                  const std::string &WeightName1,
		  const std::string &WeightName2,
		  const double Weight1Fraction,
                  const std::string &Cut = "") {
  Tree->Draw(">> elist", Cut.c_str(), "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  if(!elist) {
    std::cout << "Something's wrong!\n";
    std::cout << Cut << "\n";
    return 0.0;
  }
  Tree->SetEntryList(elist);
  double Total = 0.0, Weight1, Weight2;
  Tree->SetBranchAddress(WeightName1.c_str(), &Weight1);
  if(Weight1Fraction != 1.0) {
    Tree->SetBranchAddress(WeightName2.c_str(), &Weight2);
  }
  for(Long64_t i = 0; i < elist->GetN(); i++) {
    Tree->GetEntry(Tree->GetEntryNumber(i));
    Total += Weight1Fraction*Weight1;
    if(Weight1Fraction != 1.0) {
      Total += (1.0 - Weight1Fraction)*Weight2;
    }
  }
  Tree->SetEntryList(nullptr);
  return Total;
}

/**
 * Sum the event weights of a TTree
 * @param Tree The TTree
 * @param WeightName The name of the weights
 * @param Cut Cut to apply, or leave empty
 */
double SumWeights(TTree *Tree,
                  const std::string &WeightName = "ModelWeight",
                  const std::string &Cut = "") {
  return SumWeights(Tree, WeightName, "", 1.0, Cut);
}

/**
 * Get the number of reconstructed single tag events
 * @param TagMode The name of the tag mode
 * @param WeightName The name of the weights, or leave empty for no weighting
 */

double GetReconstructedSTEvents(
  const std::string &TagMode, const std::string &WeightName = "") {
  std::string Filename = "${BES3_ANALYSIS_PATH}/Selection/SignalMC/SingleTag/";
  Filename += TagMode + "/";
  Filename += TagMode + "_SingleTag_SignalMC.root";
  TChain Chain((TagMode + "SingleTag").c_str());
  Chain.Add(Filename.c_str());
  if(WeightName.empty()) {
    return Chain.GetEntries("abs(Run) < 50000") + 
           Chain.GetEntries("abs(Run) >= 50000")*(17.0/5.0);
  } else {
    return SumWeights(&Chain, WeightName, "abs(Run) < 50000") +
           SumWeights(&Chain, WeightName, "abs(Run) >= 50000")*(17.0/5.0);
  }
}

/**
 * Get the single tag efficiency
 * @param TagMode The name of the tag mode
 * @param WeightName The name of the weights, or leave empty for no weighting
 */
std::pair<double, double> GetSTEfficiency(
  const std::string &TagMode, const std::string &WeightName = "") {
  const double RecEvents = GetReconstructedSTEvents(TagMode, WeightName);
  const double GenEvents = 1000000.0;
  std::pair<double, double> Efficiency;
  Efficiency.first = RecEvents/GenEvents;
  Efficiency.second = Efficiency.first*(1 - Efficiency.first)/GenEvents;
  Efficiency.second = TMath::Sqrt(Efficiency.second);
  return Efficiency;
}

/**
 * Get LaTeX name of decay mode
 * Mode Name of mode
 */
std::string GetTagName(const std::string &Tag) {
  const std::map<std::string, std::string> TagNames{
    {"KKpipi", "$\\kaonp\\kaonm\\pip\\pim$"},
    {"Kpi", "$\\kaonm\\pip$"},
    {"Kpipi0", "$\\kaonm\\pip\\piz$"},
    {"Kpipipi", "$\\kaonm\\pip\\pim\\pip$"},
    {"KeNu", "$\\kaonm e^{+}\\nu_e$"},
    {"KK", "$\\kaonp\\kaonm$"},
    {"KKPartReco", "$\\kaonp\\kaonm$ part. reco."},
    {"pipi", "$\\pip\\pim$"},
    {"pipipi0", "$\\pip\\pim\\piz$"},
    {"KSpi0pi0", "$\\kshort\\piz\\piz$"},
    {"KLpi0", "$\\klong\\piz$"},
    {"KSpi0", "$\\kshort\\piz$"},
    {"KSpi0PartReco", "$\\kshort\\piz$ part. reco."},
    {"KSeta", "$\\kshort\\eta$"},
    {"KSetaPrimepipieta", "$\\kshort\\eta^\\prime_{\\pi\\pi\\eta}$"},
    {"KSetaPrimerhogamma", "$\\kshort\\eta^\\prime_{\\rho\\gamma}$"},
    {"KSpipipi0", "$\\kshort\\pip\\pim\\piz$"},
    {"KSpipi", "$\\kshort\\pip\\pim$"},
    {"KLpipi", "$\\klong\\pip\\pim$"},
    {"KSpipiPartReco", "$\\kshort\\pip\\pim$ part. reco."}
  };
  return TagNames.at(Tag);
}

/**
 * Print number with its uncertainty in LaTeX format with the correct number of decimals
 * Algorithm by Sophia Vaughan
 * Value Number to print
 * Error The uncertainty
 * Power The power to separate out
 */
std::string PrintLaTeXNumber(double Value, double Error, int Power = 0) {
    if(Power != 0) {
        Value *= TMath::Power(10, -Power);
        Error *= TMath::Power(10, -Power);
    }
    int Precision;
    if(Error >= TMath::Sqrt(10.0)) {
      Precision = 0;
    } else {
      double whole, fractional;
      fractional = std::modf(TMath::Log10(Error), &whole);
      Precision = TMath::Power(10, fractional) > TMath::Sqrt(10)/10.0 ? 1 : 2;
      Precision -= whole;
    }
    std::stringstream ss;
    ss << std::fixed << std::setprecision(Precision) << Value << " \\pm " << Error;
    std::string Line;
    std::getline(ss, Line);
    if(Power != 0) {
        Line.insert(0, "(");
        Line  += ")\\times 10^{" + std::to_string(Power) + "}";
    }
    Line.insert(0, "$");
    Line.push_back('$');
    return Line;
}

/**
 * Print a line with a variable name and a number with errors
 * VariableName Variable name (LaTeX code with $)
 * Value Number to print
 * Error The uncertainty
 * Power The power to separate out
 * BlankSpaces Number of blank spaces in front
 * Width The width of each cell
 */
void PrintLaTeXLine(const std::string &VariableName,
		    double Value,
		    double Error,
		    int Power = 0,
		    int BlankSpaces = 8,
		    int Width = 25) {
    std::string Blank = BlankSpaces == 0 ? "" : std::string(BlankSpaces, ' ');
    int width = 25;
    std::cout << Blank;
    std::cout << std::left << std::setw(Width) << VariableName << " & "
              << std::left << std::setw(Width) << PrintLaTeXNumber(Value, Error, Power) << " \\\\" << "\n";
}

/**
 * Get the bin yields in this TTree
 * Tree The TTree we're looking at
 * Reco Set to true to use the reconstructed bins
 * NumberBins The number of bins
 * WeightName If reweighted, specify the weight name
 * BinName The name of the bin variable
 */
std::map<int, double> GetBinYields(TTree *Tree,
				   bool Reco,
				   int NumberBins,
				   const std::string &WeightName = "",
				   std::string BinName = "SignalBin") {
  if(!Reco) {
    BinName += "_true";
  }
  std::map<int, double> BinYields;
  for(int Bin = -NumberBins; Bin <= NumberBins; Bin++) {
    if(Bin == 0) {
      continue;
    }
    const std::string Cut(BinName + " == " + std::to_string(Bin));
    double Events;
    if(WeightName.empty()) {
      Events = Tree->GetEntries(Cut.c_str());
    } else {
      Events = SumWeights(Tree, WeightName, Cut);
    }
    BinYields.insert({Bin, Events});
  }
  return BinYields;
}

/**
 * Get the bin yields of reconstructed signal KKpipi
 * TagMode The tag mode
 * NumberBins Number of bins
 * WeightName If reweighted, specify the weight name
 * BinName The name of the bin variable
 */
std::map<int, double> GetRecSignalBinYields(const std::string &TagMode,
					    int NumberBins,
					    const std::string &WeightName = "",
					    const std::string &BinName = "SignalBin") {
  std::string Filename = "${BES3_ANALYSIS_PATH}/Selection/SignalMC/DoubleTag/";
  Filename += TagMode + "/KKpipi_vs_" + TagMode + "_Binned_SignalMC.root";
  TChain Chain((TagMode + "DoubleTag").c_str());
  Chain.Add(Filename.c_str());
  return GetBinYields(&Chain, true, NumberBins, WeightName, BinName);
}

/**
 * Function for parsing parameters from a file into a map
 * Filename Filename, obviously
 */
std::map<std::string, double> ParseParameters(const std::string &Filename) {
  std::map<std::string, double> FittedParameters;
  std::ifstream File(Filename);
  std::string Line;
  while(std::getline(File, Line)) {
    if(Line.empty()) {
      continue;
    }
    std::stringstream ss(Line);
    double Value;
    std::string Name;
    ss >> Name >> Value;
    FittedParameters.insert({Name, Value});
  }
  File.close();
  return FittedParameters;
}

/**
 * Get the inclusive DT efficiency
 * Tag Tag mode
 * WeightName The name of the weight, otherwise leave empty
 */
std::pair<double, double> GetInclusiveDTEff(const std::string &Tag, const std::string &WeightName) {
  std::string Filename = "${BES3_ANALYSIS_PATH}/Selection/SignalMC/DoubleTag/";
  Filename += Tag + "/KKpipi_vs_" + Tag + "_Binned_SignalMC_Reweighted.root";
  TChain Chain((Tag + "DoubleTag").c_str());
  Chain.Add(Filename.c_str());
  double RecEvents, GenEvents;
  if(WeightName.empty()) {
    RecEvents = Chain.GetEntries("abs(Run) < 50000")
              + Chain.GetEntries("abs(Run) >= 50000")*(17.0/5.0);
    GenEvents = 2000000.0;
    if(Tag == "KSpipi" || Tag == "KSpipiPartReco" || Tag == "KLpipi") {
      GenEvents *= 8;
    }
  } else {
    std::string TruthFilename = "${BES3_ANALYSIS_PATH}/TruthTuples/BinnedTruthTuples/";
    TruthFilename += Tag + "/KKpipi_vs_" + Tag + "_TruthTuple_Binned_Reweighted.root";
    TChain TruthChain("TruthTuple");
    TruthChain.Add(TruthFilename.c_str());
    RecEvents = SumWeights(&Chain, WeightName, "abs(Run) < 50000")
              + SumWeights(&Chain, WeightName, "abs(Run) >= 50000")*(17.0/5.0);
    GenEvents = SumWeights(&TruthChain, WeightName, "abs(Run) < 50000")
              + SumWeights(&TruthChain, WeightName, "abs(Run) >= 50000")*(17.0/5.0);
  }
  double Eff = RecEvents/GenEvents;
  double Eff_err = TMath::Sqrt(Eff*(1.0 - Eff)/GenEvents);
  return std::make_pair(Eff, Eff_err);
}

/**
 * Calculate the covariance between two vectors
 */
template<typename T>
T Covariance(const std::vector<T> &x, const std::vector<T> &y) {
    T Mean_x = TMath::Mean(x.begin(), x.end());
    T Mean_y = TMath::Mean(y.begin(), y.end());
    T Total2 = std::inner_product(x.begin(), x.end(), y.begin(), static_cast<T>(0),
                                  std::plus<>(), [=](T a, T b) { return (a - Mean_x)*(b - Mean_y); });
    return Total2/static_cast<T>(x.size() - 1);
}

/**
 * Save covariance matrix
 * FlatCovMatrix A flattened covariance matrix
 * Filename The file to save it to
 */
void SaveCovMatrix(const std::vector<double> &FlatCovMatrix,
		   const std::string &Filename) {
  auto Size = static_cast<std::size_t>(TMath::Sqrt(FlatCovMatrix.size()));
  if(Size*Size != FlatCovMatrix.size()) {
    std::cout << "Error! Covariance matrix not square!\n";
  }
  TMatrixT<double> CovMatrix(Size, Size);
  for(std::size_t i = 0; i < Size; i++) {
    for(std::size_t j = 0; j < Size; j++) {
      CovMatrix(i, j) = FlatCovMatrix[i*Size + j];
    }
  }
  TFile File(Filename.c_str(), "RECREATE");
  File.WriteObject(&CovMatrix, "CovMatrix");
  File.Close();
}
