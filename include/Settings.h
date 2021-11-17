// Copied from GGSZ code repository
#ifndef Settings_h
#define Settings_h

#include <cstdlib>
#include <string>
#include <map>



class Settings{

public:
    Settings(
        std::string settings_name,
        std::string file_name
        );
    Settings(
        std::string settings_name="GGSZSettings_with_no_name_set"
        );
    ~Settings(){};

    std::string get_name()const {return _name;} ;

    void set_value(
        std::string key, 
        std::string val,
        std::string origin_comment="var manually set",
        bool enforce_var_already_existing=true);
    void update_from_file(
        std::string file_name,
        bool enforce_var_already_existing=true);
    void add_subsettings(
        std::string file_key,
        std::string file_val,
        bool enforce_var_already_existing = false);
    void update_subsettings_from_file(
        std::string key,
        std::string file_path,
        bool enforce_var_already_existing = false);

    // functions to get settings value 
    std::string get(
        std::string key,
        std::string default_val = "") const; // get settings value
    bool getB(std::string key) const; // get a boolean value
    double getD(std::string key) const; // get settings value in type T
    int getI(std::string key) const; // get settings value in type T

    // Overload operator to fetch subsettings
    Settings& operator [] (std::string settings_name) const;

    // check if settings contains key/subsettings
    bool contains(std::string key) const;
    bool contains_subsettings(std::string subsettings_name) const;

    // write out all settings and subsettings to file
    // (in a format that can be directly read into a new GGSZSettings instance)
    void dump_settings_to_file(
        std::string file_name,
        std::string prefix="",
        bool append_to_file=false);

private:
    void process_subsettings(
        std::string line,
        bool enforce_var_already_existing);
    void add_entry(
        std::string key,
        std::string val,
        std::string file_name,
        bool enforce_var_already_existing);
    void add_subsetting_entry(
        std::string key,
        std::string val,
        std::string file_name,
        bool enforce_var_already_existing);



    std::map<std::string, std::string> _var_map; // holds settings values
    std::map<std::string, std::string> _var_path_map; // holds settings values
    std::map<std::string, Settings *> _settings_map; // holds subsettings
    std::map<std::string, std::string> _settings_path_map; // holds subsetting paths

    std::string _name;
};


#endif
