// Martin Duy Tat 31st March 2021

#include<iostream>
#include<fstream>
#include<string>
#include<regex>
#include<stdexcept>
#include"TChain.h"
#include"RooRealVar.h"
#include"Utilities.h"
#include"InitialCuts.h"
#include"DeltaECut.h"
#include"TruthMatchingCuts.h"
#include"Settings.h"
#include"Unique.h"

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

  TCut LoadCuts(const std::string &TagMode, const std::string &TagType, bool IncludeDeltaECuts, const std::string &DataMC, const std::string &TruthMatchMode) {
    InitialCuts initialCuts(TagMode, TagType);
    TCut Cuts = initialCuts.GetInitialCuts();
    if(IncludeDeltaECuts) {
      DeltaECut deltaECut(TagMode, TagType, DataMC);
      Cuts = Cuts && deltaECut.GetDeltaECut();
    }
    if(TruthMatchMode != "") {
      TruthMatchingCuts truthMatchingCuts(TruthMatchMode, TagType);
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
    if(Mode == "Kpi" || Mode == "Kpipi0" || Mode == "Kpipipi" || Mode == "KeNu") {
      return std::unique_ptr<KKpipi_PhaseSpace>{new KKpipi_vs_Flavour_PhaseSpace(Tree, settings["BinningScheme"].getI("NumberBins"))};
    } else {
      throw std::invalid_argument(Mode + " tag mode is unknown");
    }
  }

}
