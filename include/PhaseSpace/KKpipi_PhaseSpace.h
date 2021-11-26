// Martin Duy Tat 26th November 2021
/**
 * KKpipi_PhaseSpace stores variables that are connected to a TTree and when an event is loaded it can determine the KKpipi phase space bin
 * This is assuming a \f$D^0\to KKpipi\f$ decay, but this class knows nothing of the tag side
 * This is a pure virtual class because it needs further information about the tag side to return a correct phase space bin
 */

#ifndef KKPIPI_PHASESPACE
#define KKPIPI_PHASESPACE

#include<vector>
#include<utility>
#include"TTree.h"
#include"AmplitudePhaseSpace.h"

class KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that connects the TTree branches with double tag events to the kinematic variables of the KKpipi decay
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     */
    KKpipi_PhaseSpace(TTree *Tree, int Bins);
    /**
     * Pure virtual function that returns both the signal and tag side binning
     */
    virtual std::pair<int, int> Bin() const = 0;
  protected:
    /**
     * Get the phase space bin of the \f$D^0\to KK\pi\pi\f$ decay (but obviously we know nothing about the flavour yet)
     */
    int KKpipiBin() const;
  private:
    /**
     * Vector of the kinematic variables before Kalman fit
     */
    std::vector<double> m_Momenta;
    /**
     * Vector of the kinematic variables after Kalman fit
     */
    std::vector<double> m_MomentaKalmanFit;
    /**
     * Flag that is 1 when Kalman fit was a success
     */
    int m_KalmanFitSuccess;
    /**
     * Object that calculates the phase space point of an event
     */
    AmplitudePhaseSpace m_AmplitudePhaseSpace;
};

#endif
