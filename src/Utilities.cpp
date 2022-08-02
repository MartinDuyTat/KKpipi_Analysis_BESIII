// Martin Duy Tat 31st March 2021

#include<iostream>
#include<fstream>
#include<string>
#include<regex>
#include<stdexcept>
#include"TChain.h"
#include"TTree.h"
#include"TEntryList.h"
#include"RooRealVar.h"
#include"Utilities.h"
#include"InitialCuts.h"
#include"DeltaECut.h"
#include"TruthMatchingCuts.h"
#include"Settings.h"
#include"Unique.h"
#include"PhaseSpace/KKpipi_PhaseSpace.h"
#include"PhaseSpace/KKpipi_vs_CP_PhaseSpace.h"
#include"PhaseSpace/KKpipi_vs_Flavour_PhaseSpace.h"
#include"PhaseSpace/KKpipi_vs_K0hh_PhaseSpace.h"

namespace Utilities {
  void LoadChain(TChain *Chain, int NumberFiles, const std::string &Filename, const std::string &TreeName) {
    std::cout << "Initializing TChain with files...\n";
    if(TreeName != "") {
      Chain->SetName(TreeName.c_str());
    }
    if(NumberFiles == 0) {
      std::cout << "Need more than 0 input files...\n";
      return;
    }
    for(int i = 0; i < NumberFiles; i++) {
      Chain->Add((Filename + std::to_string(i) + ".root").c_str());
    }
    std::cout << "ROOT files added to TChain\n";
    return;
  }

  void LoadChain(TChain *Chain, const std::string &Filename, const std::string &TreeName) {
    std::cout << "Initializing TChain with files...\n";
    if(TreeName != "") {
      Chain->SetName(TreeName.c_str());
    }
    std::ifstream Infile(Filename);
    std::string line;
    while(std::getline(Infile, line)) {
      Chain->Add(line.c_str());
    }
    Infile.close();
    std::cout << "ROOT files added to TChain\n";
    return;
  }

  TCut LoadCuts(const std::string &SignalMode,
		const std::string &TagMode,
		const std::string &TagType,
		const std::string &DataMC,
		bool IncludeDeltaECuts,
		bool TruthMatch,
		bool KKpipiPartReco) {
    std::string Mode = TagMode;
    if(Mode.find("_to_") != std::string::npos) {
      Mode = Mode.substr(Mode.find("_to_") + 4);
    }
    InitialCuts initialCuts(Mode, TagType, KKpipiPartReco);
    TCut Cuts = initialCuts.GetInitialCuts();
    if(IncludeDeltaECuts) {
      DeltaECut deltaECut(Mode, TagType, DataMC, KKpipiPartReco);
      Cuts = Cuts && deltaECut.GetDeltaECut();
    }
    if(TruthMatch) {
      TruthMatchingCuts truthMatchingCuts(SignalMode, TagMode);
      Cuts = Cuts && truthMatchingCuts.GetTruthMatchingCuts();
    }
    return Cuts;
  }

  void replace_env_variables(std::string& text) {
    static std::regex env( "\\$\\{([^}]+)\\}" );
    std::smatch match;
    while ( std::regex_search( text, match, env ) ) {
        const char * s = getenv( match[1].str().c_str() );
        const std::string var( s == NULL ? "" : s );
        text.replace( match.position(0), match.length(0), var );
    }
  }

  Settings parse_args(int argc, char** argv) {
    // Parses options and returns a GGSZSettings instance with the corresponding options
    if (argc < 2) {
      exit(-1);
    }
    std::string default_settings_file = argv[1];
    Settings settings("general_settings", default_settings_file);
    for (int i = 1; i < argc; i++){
      // First a loop simply detecting help calls
      std::string arg = argv[i];
      if (arg=="-h") {
	exit(0);
      }
    }
    // If no help calls were made, process arguments
    for (int i = 1; i < argc; i++){
      std::string arg = argv[i];
      if (arg=="-es"){
	std::string extra_settings_file = argv[i+1];
	settings.update_from_file(extra_settings_file);
	i = i + 1; // don't process arg that is settings file
      }
      if (arg=="-o"){
	std::string key = argv[i+1];
	std::string val = argv[i+2];
	bool enforce_key_already_existing = true;
	settings.set_value(key, val, "cmd line", enforce_key_already_existing);
	i = i + 2; // don't process args with key/value
      }
      
      if (arg=="-f"){
	std::string key = argv[i+1];
	std::string val = argv[i+2];
	bool enforce_key_already_existing = false; // can be nice to include NEW option keys for syst
	settings.update_subsettings_from_file(key, val, enforce_key_already_existing);
	i = i + 2; // don't process args with key/value
      }   
    }
    return settings;
  }

  RooRealVar* load_param(const Settings& settings, const std::string& name) {
    double init_val = settings.getD(name); // throws error if var not included
    double low = init_val;
    double high = init_val;
    if (settings.contains(name+"_limL"))
        low = settings.getD(name+"_limL");
    if (settings.contains(name+"_limU"))
        high = settings.getD(name+"_limU");
    RooRealVar * v = Unique::create<RooRealVar*>(name, "", init_val, low, high);
    // Handle limits: variables without limL (limU) are set to have no lower (upper) limit
    if (!settings.contains(name+"_limL"))
        v->removeMin();
    if (!settings.contains(name+"_limU"))
        v->removeMax();
    // If limits are equal, set constant
    if (low==high){
        v->setConstant();
    }
    // If <var_name>_float is in settings, it determines whether param. is floated
    if (settings.contains(name+"_float")){
        v->setConstant(!settings.getB(name+"_float"));
    }
    return v;
  }

  std::vector<std::string> ConvertStringToVector(std::string List) {
    std::replace(List.begin(), List.end(), ',', ' ');
    std::stringstream ss(List);
    std::vector<std::string> ListOfElements;
    for(std::string i; ss >> i;) {
      ListOfElements.push_back(i);
    }
    return ListOfElements;
  }

  std::string ReplaceString(const std::string &String, const std::string &From, const std::string &To) {
    return std::regex_replace(String, std::regex(From), To);
  }

  std::unique_ptr<KKpipi_PhaseSpace> GetPhaseSpaceBinning(const Settings &settings, TTree *Tree) {
    std::string Mode = settings.get("Mode");
    int NumberBins = settings["BinningScheme"].getI("NumberBins");
    bool Bin_rec = settings.getB("Bin_reconstructed");
    bool Bin_truth = settings.getB("Bin_truth");
    bool KSKKBkg = settings.contains("KSKKBackground") && settings.getB("KSKKBackground");
    if(GetTagType(Mode) == "Flavour") {
      return std::unique_ptr<KKpipi_PhaseSpace>{new KKpipi_vs_Flavour_PhaseSpace(Tree, NumberBins, Bin_rec, Bin_truth, KSKKBkg)};
    } else if(GetTagType(Mode) == "SCMB") {
      bool KStoKLBackground = Mode.substr(0, 2) == "KL" && settings.get("TagMode").substr(0, 2) == "KS";
      return std::unique_ptr<KKpipi_PhaseSpace>{new KKpipi_vs_K0hh_PhaseSpace(Tree, NumberBins, Bin_rec, Bin_truth, Mode, KSKKBkg, settings.contains("KKpipiPartReco") && settings.getB("KKpipiPartReco"), KStoKLBackground)};
    } else if(GetTagType(Mode) == "CP") {
      return std::unique_ptr<KKpipi_PhaseSpace>{new KKpipi_vs_CP_PhaseSpace(Tree, NumberBins, Bin_rec, Bin_truth, KSKKBkg)};
    } else {
      return nullptr;
    }
  }

  std::string GetTagType(const std::string &Mode) {
    if(Mode == "Kpi" || Mode == "Kpipi0" || Mode == "Kpipipi" || Mode == "KeNu") {
      return "Flavour";
    } else if(Mode == "KSpipi" || Mode == "KSKK" || Mode == "KLpipi" || Mode == "KLKK") {
      return "SCMB";
    } else {
      const std::vector<std::string> CPModes{"KK", "pipi", "KSpi0pi0", "pipipi0", "KLpi0", "KSpi0", "KSeta", "KSetaPrimepipieta", "KSetaPrimerhogamma", "KSpipipi0", "KLpi0pi0"};
      if(std::find(CPModes.begin(), CPModes.end(), Mode) != CPModes.end()) {
	return "CP";
      } else {
	throw std::invalid_argument(Mode + " tag mode is unknown");
      }
    }
  }

  double SumWeights(TTree *Tree, const std::string &WeightName, const std::string &Cut) {
    TTree *ClonedTree = Tree->CloneTree();
    ClonedTree->Draw(">> elist", Cut.c_str(), "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    ClonedTree->SetEntryList(elist);
    double Total = 0.0, Weight;
    ClonedTree->SetBranchAddress(WeightName.c_str(), &Weight);
    for(Long64_t i = 0; i < elist->GetN(); i++) {
      ClonedTree->GetEntry(ClonedTree->GetEntryNumber(i));
      Total += Weight;
    }
    return Total;
  }

  std::string GetTagNameLaTeX(const std::string &Tag) {
    if(Tag == "KKpipi") {
      return "K^{+}K^{#minus}#pi^{+}#pi^{#minus}";
    } else if(Tag == "KK") {
      return "K^{+}K^{#minus}";
    } else if(Tag == "pipi") {
      return "#pi^{+}#pi^{#minus}";
    } else if(Tag == "pipipi0") {
      return "#pi^{+}#pi^{#minus}#pi^{0}";
    } else if(Tag == "KSpi0pi0") {
      return "K_{S}^{0}#pi^{0}#pi^{0}";
    } else if(Tag == "KLpi0") {
      return "K_{L}^{0}#pi^{0}";
    } else if(Tag == "KSpi0") {
      return "K_{S}^{0}#pi^{0}";
    } else if(Tag == "KSeta") {
      return "K_{S}^{0}#eta";
    } else if(Tag == "KSetaPrimepipieta") {
      return "K_{S}^{0}#eta'(#pi^{+}#pi^{#minus}#eta)";
    } else if(Tag == "KSetaPrimerhogamma") {
      return "K_{S}^{0}#eta'(#pi^{+}#pi^{-}#gamma)";
    } else if(Tag == "KSpipipi0") {
      return "K_{S}^{0}#omega";
    } else if(Tag == "KSpipi" || Tag == "KSpipiPartReco") {
      return "K_{S}^{0}#pi^{+}#pi^{#minus}";
    } else if(Tag == "KLpipi") {
      return "K_{L}^{0}#pi^{+}#pi^{#minus}";
    } else {
      throw std::invalid_argument("Unknown tag: " + Tag);
    }
  }

}
