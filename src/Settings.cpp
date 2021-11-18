#include "Settings.h"
#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <set>
#include "boost/algorithm/string.hpp"



Settings::Settings(
    std::string settings_name,
    std::string file_name
    ){
    // Constructor that loads settings from a file
    _name = settings_name;
    this->update_from_file(file_name, false);
}


Settings::Settings(
    std::string settings_name){
    // Constructor that loads no settings (yet)
    _name = settings_name;
}



void Settings::update_from_file(
    std::string file_name,
    bool enforce_var_already_existing
    ){
    // Load settings from a file
    // if a settings value already exists, it is overwritten (and an INFO statement printed)
    // if no value exists an error is thrown if "enforce_already_existing = true" to catch misspelled settings overrides
    // if file holds same settings more than once, an error is thrown (as this is likely a bug)
    std::ifstream file;
    std::string line;
    file.open(file_name.c_str());

    if (!file.is_open()){
        std::cerr << "Error opening file: " << file_name << " in settings: " << _name << "\n";
        throw std::runtime_error("Error opening settings file in settings!");
    }

    std::set<std::string> added_keys;

    while (file.good() && !file.eof()){
        std::getline(file, line);
        line = line.substr(0, line.find("*")); // ignore comments
        boost::trim(line); // remove any spaces
        if (line.size()==0) continue;
        int space_pos = line.find(' '); // find break
        if (space_pos == (int)std::string::npos){ // line must have key and val, but hasn't
            std::cout << "Skipping '" << line << "' in settings: " << _name << ". No value given\n";
            continue;
        }
        std::string key = line.substr(0, space_pos);
        std::string val = line.substr(1 + space_pos, line.length()-space_pos-1);

        if (key=="file"){
            // for special key 'file' load subsettings
            this->process_subsettings(val, enforce_var_already_existing);
            continue;
        }

        if (added_keys.find(key)!=added_keys.end()){
            std::cerr << "key '" << key << "' found twice in file '" << file_name << "'!\n";
            throw std::runtime_error("same key in settings file twice!");
        }
        added_keys.insert(key);

        this->add_entry(key, val, file_name, enforce_var_already_existing);


    }

    file.close();
}

void Settings::add_entry(
    std::string key,
    std::string val,
    std::string file_name,
    bool enforce_var_already_existing){

    if (key.find("/") != std::string::npos){
        this->add_subsetting_entry(key, val, file_name, enforce_var_already_existing);
        return;
    }

    if (enforce_var_already_existing && !this->contains(key)){
        std::cerr << "key '" << key << "' does not exist in settings '" << _name << "'! CANNOT UPDATE!" << " ('" << file_name << "')\n";
        throw std::runtime_error("Trying to update non-existing variable!");
    }

    if (this->contains(key)){
        std::cout << "overwriting '" << key << "' in settings '" << _name << "'" << "        " << this->get(key) << " -> " << val << " ('" << file_name << "')\n";
    } 

    Utilities::replace_env_variables(val);

    _var_map[key] = val;
    _var_path_map[key] = file_name;

}

void Settings::add_subsetting_entry(
    std::string key,
    std::string val,
    std::string file_name,
    bool enforce_var_already_existing){

    int slash_pos = key.find('/');
    std::string settings_name = key.substr(0, slash_pos);
    std::string subsettings_key = key.substr(1 + slash_pos, key.length()-slash_pos-1);
    if (!this->contains_subsettings(settings_name)){
        std::string subsettings_name =  settings_name;
        std::cout << "adding new subsettings '" << subsettings_name << "' to settings '" << _name << "\n";
        _settings_map[settings_name] = new Settings(subsettings_name) ;
        _settings_path_map[settings_name] = file_name;
    }
    _settings_map[settings_name]->add_entry(subsettings_key, val, file_name, enforce_var_already_existing);
}

void Settings::process_subsettings(
    std::string val,
    bool enforce_var_already_existing){

    int file_space_pos = val.find(' ');
    if (file_space_pos == (int)std::string::npos){ // line must have key and val, but hasn't
        std::cout << "Skipping file '" << val << "' in settings: " << _name << ". No file path given\n";
        return;
    }
    std::string file_key = val.substr(0, file_space_pos);
    std::string file_val = val.substr(1 + file_space_pos, val.length()-file_space_pos-1);
    

    add_subsettings(file_key, file_val, enforce_var_already_existing);

}

void Settings::add_subsettings(
    std::string file_key,
    std::string file_val,
    bool enforce_var_already_existing){

    Utilities::replace_env_variables(file_val);
    std::string subsettings_name = file_key;
    if (!this->contains_subsettings(file_key)){
        // subsettings does not exist, so create it
        std::cout << "Added subsettings '" << file_key << "' from file '" << file_val << "' to settings: " << _name << "\n";
        _settings_map[file_key] = new Settings(subsettings_name, file_val) ;
        _settings_path_map[file_key] = file_val;
    } else {
        // subsettings does exist, so update it
        std::cout << "UPDATING subsettings '" << file_key << "' from file '" << file_val << "' in settings: " << _name << "\n";
        _settings_map[file_key]->update_from_file(file_val, enforce_var_already_existing);
        _settings_path_map[file_key] = _settings_path_map[file_key] + ", " + file_val;
    }
}

void Settings::update_subsettings_from_file(
    std::string key,
    std::string file_path,
    bool enforce_already_existing){

    if (key.find("/") != std::string::npos){
        // Trying to update a file further down the tree: simply pass down!
        int slash_pos = key.find('/');
        std::string settings_name = key.substr(0, slash_pos);
        std::string subsettings_name = key.substr(1 + slash_pos, key.length()-slash_pos-1);
        _settings_map[settings_name]->update_subsettings_from_file(subsettings_name, file_path, enforce_already_existing);
        return;
    }

    if (!this->contains_subsettings(key)){
        std::cerr << "Cannot update non-existing subsettings!\n";
        throw std::runtime_error("Trying to update non-existent subsettings: '" + key + "' in settings: " + _name);
    }

    _settings_map[key]->update_from_file(file_path, enforce_already_existing);
    return;

}

void Settings::set_value(
    std::string key,
    std::string value,
    std::string origin_comment,
    bool enforce_var_already_existing){

    this->add_entry(key, value, origin_comment, enforce_var_already_existing);
 // origin of var value stored for easy reproduceablity
}

bool Settings::contains(
    std::string key) const{
    return _var_map.find(key)!=_var_map.end();
}

bool Settings::contains_subsettings(
    std::string key) const{
    return _settings_map.find(key)!=_settings_map.end();
}

Settings& Settings::operator[] (
    std::string settings_name) const{

    if (!this->contains_subsettings(settings_name)){
        std::cerr << "Trying to access non-existing subsettings '" << settings_name << "' in settings: " << _name << "\n";
        throw std::runtime_error("Trying to access non-existing subsettings in settings");
    }
    return *(_settings_map.at(settings_name));
}

std::string Settings::get(
    std::string key,
    std::string default_val) const {

    if (contains(key)) return _var_map.at(key);
    if (default_val != "") return default_val;
    // throw error if key does not exist and no default value supplied
    std::cerr << "Trying to access non-existing key '" << key << "' in settings: " << _name << "\n";
    throw std::runtime_error("Trying to access non-existing key in settings");
}

bool Settings::getB(
    std::string key) const{

    if (this->contains(key)){
        std::string val = this->get(key);
        boost::algorithm::to_lower(val);
        return (val == "1" || val=="true");
    }
    // throw error if key does not exist
    std::cerr << "Trying to access non-existing key '" << key << "' in as BOOL settings: " << _name << "\n";
    throw std::runtime_error("Trying to access non-existing key as BOOL in settings: " + _name);

}

double Settings::getD(
    std::string key) const {

    if (this->contains(key)){
        return atof(this->get(key).c_str());
    }
    std::cerr << "Trying to access non-existing key '" << key << "' in as T settings: " << _name << "\n";
    throw std::runtime_error("Trying to access non-existing key as T in settings: " + _name);
}

int Settings::getI(
    std::string key) const {

    if (this->contains(key)){
        return atoi(this->get(key).c_str());
    }
    std::cerr << "Trying to access non-existing key '" << key << "' in as T settings: " << _name << "\n";
    throw std::runtime_error("Trying to access non-existing key as T in settings: " + _name);
}

void Settings::dump_settings_to_file(
    std::string file_name,
    std::string prefix,
    bool append_to_file){

    std::ofstream ofile;
    if (append_to_file) ofile.open(file_name.c_str(), std::ios::app);
    else ofile.open(file_name.c_str());


    // First write all variables in this settings
    ofile << "* Settings in : " << _name << std::endl;
    for (auto const& v : _var_map){
        ofile << prefix << v.first << " " << v.second << " * from: " << _var_path_map[v.first] << std::endl;
    }
    ofile << std::endl;
    ofile.close();


    // Then make all subsettings do the same in turn
    // now using a prefix when writing, so that the whole file "file_name"
    // can be read and original settings structure recovered
    for (auto const& s : _settings_map){
        std::string new_prefix = prefix + s.second->get_name() + "/";
        s.second->dump_settings_to_file(file_name, new_prefix, true);
    }

}

