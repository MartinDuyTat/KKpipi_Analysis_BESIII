// Martin Duy Tat 27th April 2021
/**
 * CutsFromFile is a class that loads cuts from a .cut file inside the library
 * These files are either in the directory InitalCuts or TruthMatchingCuts and are mode-dependent
 * The variables in the cuts text file must have the format "SignalTagVariable", and the "SignalTag" part will be replaced appropriately
 */

#ifndef CUTSFROMFILE
#define CUTSFROMFILE

#include<string>
#include"TCut.h"

class CutsFromFile {
  public:
    /**
     * Helper function that replaces the "SignalTag" part of variables with the correct tag type
     * @param Input string with cuts
     * @param TagSide Either "Signal" (KKpipi) or "Tag" (any other tag mode) for double tags, leave blank for single tags
     * @return Returns a string with the correct tag mode
     */
    std::string ReplaceTagMode(std::string Cuts, const std::string &TagSide) const;
    /**
     * Helper function that reads the cuts from a file
     * @param TagMode "KKpipi", "Kpi", etc
     * @param TagSide Either "Signal" or "Tag" for double tags, leave blank for single tags
     * @return Returns a TCut object with the cuts for that tag side
     */
    TCut GetTagCuts(const std::string &TagMode, const std::string &TagSide) const;
    /**
     * Constructor that takes in the filename of the cuts and produces the correct TCut object
     * @param Filename Filename of file with the cuts
     * @param TagMode "KKpipi", "Kpi", etc
     * @param TagSide Either "Signal" or "Tag" for double tags, leave blank for single tags
     */
    CutsFromFile(const std::string &Filename, const std::string &TagSide);
    /**p
     * Function that returns the cuts
     */
    TCut GetCuts() const;
  private:
    /**
     * The cuts read from file
     */
    TCut m_Cuts;
};

#endif
