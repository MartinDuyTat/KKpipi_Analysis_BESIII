// Martin Duy Tat 31st March 2021
/**
 * InitialCuts is a class that loads the mode-dependent initial cuts from text files in the directory InitialCuts and puts them together
 * The variables in the cuts text file must have the format "SignalTagVariable", and the "SignalTag" part will be replaced appropriately
 */

#ifndef INITIALCUTS
#define INITIALCUTS

#include<string>
#include"TCut.h"

class InitialCuts {
  public:
    /**
     * Constructor that sets up all the variables
     * @param TagMode The name of the tag mode we want to apply the cuts to
     * @param TagType Leave blank for single tags, input "Signal" or "Tag" for the signal or tag side of a double tag
     */
    InitialCuts(const std::string &TagMode, const std::string &TagType = std::string());
    /**
     * Function that returns the complete initial cut
     */
    TCut GetInitialCuts() const;
  private:
    /**
     * The tag mode of interest
     */
    std::string m_TagMode;
    /**
     * The tag type of interest, "ST" for single tags, "DT" for double tag
     */
    std::string m_TagType;
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
    TCut GetTagCuts(const std::string &TagMode, const std::string &TagSide = std::string()) const;
};

#endif
