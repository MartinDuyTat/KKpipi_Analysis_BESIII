// Martin Duy Tat 1st December 2021
/**
 * DCS_Parameters is a class that contains the \f$r_D\f$ ratio between favoured and DCS D decays, the coherence factor \f$R\f$ and the strong phase difference \f$\delta_D\f$ between the favoured and DCS decay
 * Their uncertainties and correlations are accounted for
 */

#ifndef DCS_PARAMETERS
#define DCS_PARAMETERS

#include<vector>
#include"Settings.h"
#include"uncertainties/ureal.hpp"

class DCS_Parameters {
  public:
    /**
     * Constructor that reads the parameters and their uncertainties and correlations from the settings
     * @param settings The settings object
     */
    DCS_Parameters(const Settings &settings);
    /**
     * Get the DCS parameters, in the order \f$r_D\f$, \f$R\f$, \f$\delta_D\f$
     */
    std::vector<uncertainties::udouble> GetDCSParameters() const;
  private:
    /**
     * The mean value of the DCS parameters
     */
    std::vector<double> m_DCS;
    /**
     * The covariance matrix of the DCS parameters, in a row-major format
     */
    std::vector<double> m_DCS_Cov;
};

#endif
