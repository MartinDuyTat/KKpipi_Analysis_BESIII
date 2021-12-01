// Martin Duy Tat 1st December 2021
/**
 * Ki is a class that stores the \f$K_i\f$ parameters
 */

#ifndef KI
#define KI

#include<vector>
#include"Settings.h"

class Ki {
  public:
    /**
     * Constructor that loads the \f$K_i\f$ from the settings
     * @param settings The settings for \f$K_i\f$
     */
    Ki(const Settings &settings);
    /**
     * Get \f$K_i\f$
     * @param Bin Bin number (i)
     */
    double Get_Ki(int Bin) const;
  private:
    /**
     * Vector of \f$K_i\f$ for positive i
     */
    std::vector<double> m_Ki;
    /**
     * Vector of \f$K_i\f$ for negative i
     */
    std::vector<double> m_Kbari;
};

#endif
