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
     * @param TagType "ST" for single tag and "DT" for double tag
     */
    TruthMatchingCuts(const std::string &TagType);
    /**
     * Function for obtaining the truth maching requirement cut
     */
    TCut GetTruthMatchingCuts() const;
  private:
    /**
     * "ST" for single tag and "DT" for double tag
     */
    std::string m_TagType;
};

#endif
