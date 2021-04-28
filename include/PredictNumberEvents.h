// Martin Duy Tat 28th April 2021
/**
 * PredictNumberEvents is a class for predicting the number of events in each bin, given the strong phase parameters and the single tag yields
 * The () operator is overloaded so that this class can be used as a functor
 */

#ifndef PREDICTNUMBEREVENTS
#define PREDICTNUMBEREVENTS

#include<vector>
#include"DDecayParameters.h"

class PredictNumberEvents {
  public:
    /**
     * Constructor that takes in a DDecayParameters object with the hadronic parameters, single tag yields and CP of the tag
     * @param DDParameters DDecayParameters object with all the \f$D\f$ hadronic parameters ready
     * @param DTFlavourYield Double tag yield of the flavour tags
     * @param STFlavourYield Single tag yield of the flavour tags (not required for flavour yields)
     * @param STYieldTag Single tag yield of the tag (not required for flavour yields)
     * @param CP CP value of the tag, \f$\pm 1\f$, or 0 for a flavour tag
     */
    PredictNumberEvents(const DDecayParameters &DDParameters, double DTFlavourYield, double STFlavourYield = 0.0, double STTagYield = 0.0, int CP = 0);
    /**
     * Function that returns the expected double tagged yield in bin \f$i\f$ for \f$CP = +1\f$ tags
     * @param i Bin number \f$i\f$
     */
    double GetCPEvenTagYield(int i) const;
    /**
     * Function that returns the expected double tagged yield in bin \f$i\f$ for \f$CP = -1\f$ tags
     * @param i Bin number \f$i\f$
     */
    double GetCPOddTagYield(int i) const;
    /**
     * Function that returns the expected double tagged yield in bin \f$i\f$ for flavour tags
     * @param i Bin number \f$i\f$
     */
    double GetFlavourTagYield(int i) const;
    /**
     * Overload of () operator so that this class can be used as a functor to obtain the double tagged yields as a function of bin number
     * @param i Bin number \f$i\f$
     */
    double operator()(int i) const;
    /**
     * Function that returns the number of bins in the binning scheme
     */
    int GetNBins() const;
  private:
    /**
     * The fractional yield of \f$D^0\f$ in bin \f$i\f$
     */
    std::vector<double> m_K;
    /**
     * The fractional yield of \f$\bar{D^0}\f$ in bin \f$i\f$
     */
    std::vector<double> m_Kbar;
    /**
     * The amplitude averaged cosine of the strong phase in bin \f$i\f$
     */
    std::vector<double> m_ci;
    /**
     * The amplitude averaged sine of the strong phase in bin \f$i\f$
     */
    std::vector<double> m_si;
    /**
     * Number of bins in the binning scheme
     */
    int m_NBins;
    /**
     * The double tag yield of the flavour tag modes
     */
    double m_DTFlavourYield;
    /**
     * The single tag yield of the flavour tag modes
     */
    double m_STFlavourYield;
    /**
     * The single tag yield of the tag mode
     */
    double m_STTagYield;
    /**
     * The \f$CP\f$ of the tag mode
     */
    int m_CP;
};

#endif
