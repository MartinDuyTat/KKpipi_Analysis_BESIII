// Martin Duy Tat 31st March 2021
/**
 * InitialCuts is a class that loads the mode-dependent initial cuts from text files in the directory InitialCuts and puts them together
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
     * The tag type of interest, empty for single tags, "Signal" or "Tag" for signal or tag side of a double tag
     */
    std::string m_TagType;
};

#endif
