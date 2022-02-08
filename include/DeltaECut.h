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
     * @param DataMC "MC" or "Data" (cuts are different for MC and data)
     * @param KKpipiPartReco Set to true for partially reconstructed KKpipi
     */
    DeltaECut(const std::string &TagMode, const std::string &TagType, const std::string &DataMC, bool KKpipiPartReco = false);
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
     * String that is either "MC" or "Data"
     */
    std::string m_DataMC;
    /**
     * Flag that is true if the KKpipi signal tag is partially reconstructed
     */
    bool m_KKpipiPartReco;
    /**
     * Helper function for reading the cuts from a file
     * @param TagMode Tag mode of interest
     * @param TagSide "Signal" or "Tag" for double tags, blank for single tags
     */
    TCut GetDeltaECutFromFile(const std::string &TagMode, const std::string &TagSide = std::string()) const;
    
};

#endif
