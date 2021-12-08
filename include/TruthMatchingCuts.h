// Martin Duy Tat 4th April 2021
/**
 * TruthMatchingCuts is a class that creates the truth matching cuts, meaning requiring that all (charged) particles have been assigned the correct particle ID and the originate from the same \f$D\f$ meson
 */

#ifndef TRUTHMATCHINGCUTS
#define TRUTHMATCHINGCUTS

#include<string>
#include"TCut.h"

class TruthMatchingCuts {
  public:
    /**
     * Constructor that sets the tag type
     * @param SignalMode "KKpipi", "KK", "pipi", etc, or "KSKK_to_KKpipi" for background modes
     * @param TagMode "KKpipi", "KK", "pipi", etc, or "Kpi_to_KK" for background modes, leave blank for single tags
     */
    TruthMatchingCuts(const std::string &SignalMode, const std::string &TagMode = "");
    /**
     * Function for obtaining the truth maching requirement cut from a file
     * If file doesn't exist, some standard truth matching is applied
     */
    TCut GetTruthMatchingCuts() const;
  private:
    /**
     * Signal mode ("KKpipi", "KK", "pipi", etc)
     */
    std::string m_SignalMode;
    /**
     * Tag mode ("KKpipi", "KK", "pipi", etc)
     */
    std::string m_TagMode;
    /**
     * Helper function that checks if truth matching cut is empty, and if it's empty a standard cut is put in
     * @param Cut Cut we want to check
     * @param TagSide "Signal" or "Tag"
     */
    TCut CheckEmptyCut(TCut Cut, const std::string &TagSide = "") const;
};

#endif
