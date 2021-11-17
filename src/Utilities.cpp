// Martin Duy Tat 31st March 2021

#include<iostream>
#include<fstream>
#include<string>
#include<regex>
#include"TChain.h"
#include"Utilities.h"
#include"InitialCuts.h"
#include"DeltaECut.h"
#include"TruthMatchingCuts.h"
#include"Settings.h"

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

  TCut LoadCuts(const std::string &CutType, const std::string &TagMode, const std::string &TagType, const std::string &DataMC) {
    InitialCuts initialCuts(TagMode, TagType);
    if(CutType == "DeltaECuts") {
      DeltaECut deltaECut(TagMode, TagType, DataMC);
      return initialCuts.GetInitialCuts() && deltaECut.GetDeltaECut();
    } else if(CutType == "NoDeltaECuts") {
      return initialCuts.GetInitialCuts();
    } else if(CutType == "TruthMatchingCuts") {
      TruthMatchingCuts truthMatchingCuts(TagMode, TagType);
      DeltaECut deltaECut(TagMode, TagType, DataMC);
      return truthMatchingCuts.GetTruthMatchingCuts() && deltaECut.GetDeltaECut() && initialCuts.GetInitialCuts();
    } else {
      return TCut();
    }
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


}
