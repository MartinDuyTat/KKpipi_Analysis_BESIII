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
#include"GeneratorKinematics.h"

class KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that connects the TTree branches with double tag events to the kinematic variables of the KKpipi decay
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     * @param ReconstructedBins Set to true to calculate the reconstructed bins
     * @param TrueBins Set to true to calculate the true bins
     */
    KKpipi_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins);
    /**
     * Pure virtual function that returns both the signal and tag side binning
     */
    virtual std::pair<int, int> Bin() const = 0;
    /**
     * Pure virtual function that returns both the signal and tag side true binning
     */
    virtual std::pair<int, int> TrueBin() = 0;
    /**
     * Get the difference between the reconstructed and true four-momenta of the D daughters
     * One must call TrueBin() first before calling this function!
     */
    std::vector<double> GetMomentumResolution() const;
  protected:
    /**
     * Get the phase space bin of the \f$D^0\to KK\pi\pi\f$ decay (but obviously we know nothing about the flavour yet)
     */
    int KKpipiBin() const;
    /**
     * Get the true phase space bin of the \f$D^0\to KK\pi\pi\f$ decay (but obviously we know nothing about the flavour yet)
     */
    int TrueKKpipiBin();
    /**
     * Struct containing the generator kinematics used to determine the true bin
     */
    GeneratorKinematics m_TrueKinematics;
    /**
     * Function that saves the index of the signal and tag D mesons in the list of particle IDs in m_TrueKinematics
     */
    void FindDIndex();
  private:
    /**
     * Flag that is 1 when Kalman fit was a success
     */
    int m_KalmanFitSuccess;
    /**
     * Vector of the kinematic variables before Kalman fit
     */
    std::vector<double> m_Momenta;
    /**
     * Vector of the kinematic variables after Kalman fit
     */
    std::vector<double> m_MomentaKalmanFit;
    /**
     * Vector of the true variables
     */
    std::vector<double> m_TrueMomenta;
    /**
     * Object that calculates the phase space point of an event
     */
    AmplitudePhaseSpace m_AmplitudePhaseSpace;
    /**
     * Set the branch addresses of reconstructed variables
     */
    void SetBranchAddresses_Rec(TTree *Tree);
    /**
     * Set the branch addresses of truth variables
     */
    void SetBranchAddresses_True(TTree *Tree);
};

#endif
