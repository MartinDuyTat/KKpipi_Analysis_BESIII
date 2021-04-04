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
     * @param TagType "ST" for single tags, "DT" double tags
     */
    DeltaECut(const std::string &TagMode, const std::string &TagType);
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
     * The tag type of interest, "ST" for single tags, "DT" for double tags
     */
    std::string m_TagType;
    /**
     * Helper function for reading the cuts from a file
     * @param TagMode Tag mode of interest
     * @param TagSide "Signal" or "Tag" for double tags, blank for single tags
     */
    TCut GetDeltaECutFromFile(const std::string &TagMode, const std::string &TagSide = std::string()) const;
    
};

#endif
