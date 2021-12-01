// Martin Duy Tat 1st December 2021
/**
 * cisi is a class that stores the \f$c_i\f$ and \f$s_i\f$ parameters
 */

#ifndef CISI
#define CISI

#include<vector>
#include"Settings.h"

class cisi {
  public:
    /**
     * Constructor that loads the \f$c_i\f$ and \f$s_i\f$ from the settings
     * @param settings The settings for \f$c_i\f$ and \f$s_i\f$
     */
    cisi(const Settings &settings);
    /**
     * Get \f$c_i\f$
     * @param Bin Bin number (i)
     */
    double Get_ci(int Bin) const;
    /**
     * Get \f$s_i\f$
     * @param Bin Bin number (i)
     */
    double Get_si(int Bin) const;
  private:
    /**
     * Vector of \f$c_i\f$
     */
    std::vector<double> m_ci;
    /**
     * Vector of \f$s_i\f$
     */
    std::vector<double> m_si;
};

#endif
