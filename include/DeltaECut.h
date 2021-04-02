// Martin Duy Tat 2nd April 2021
/**
 * DeltaECut is a class that loads the mode-dependent \f$\Delta E\f$ cuts from text files in the directory DeltaECuts and puts them together
 * The text files must have two space-separated numbers, being the lower and upper \f$\Delta E\f$ cut
 */

#ifndef DELTAECUT
#define DELTAECUT

#include<string>
#include"TCut.h"

class DeltaECut {
  public:
    /**
     * Constructor that sets up all the variables
     * @param TagMode The name of the tag mode we want to apply the cuts to
     * @param TagType Leave blank for single tags, input "Signal" or "Tag" for the signal or tag side of a double tag
     */
    DeltaECut(const std::string &TagMode, const std::string &TagType = std::string());
    /**
     * Function that returns the complete initial cut
     */
    TCut GetDeltaECut() const;
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
