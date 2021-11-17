// Martin Duy Tat 28th April 2021
/**
 * DoubleTagYield is a class for determining the double tagged yield using sideband subtraction
 */

#ifndef DOUBLETAGYIELD
#define DOUBLETAGYIELD

#include<vector>
#include"TTree.h"
#include"TCut.h"
#include"AmplitudePhaseSpace.h"

class DoubleTagYield {
  public:
    /**
     * Constructor that sets up the binning scheme and takes in the TTree with double tagged events
     * @param NBins Number of bins in binning scheme
     * @param Cuts Cuts to apply before counting events
     */
    DoubleTagYield(int NBins, TCut Cuts = TCut());
    /**
     * Helper function that determines which region of the beam constrained mass space the event belongs to
     * @param SigmalMBC Beam constrained mass of signal side (\f$KK\pi\pi\f$)
     * @param TagMBC Beam constrained mass of tag side (\f$KK\f$, \f$\pi\pi\f$, \f$K\pi\f$, etc)
     */
    char DetermineMBCRegion(double SignalMBC, double TagMBC) const;
    /**
     * Helper function that sets the correct branch addresses for momenta of daughter particles
     */
    void SetDaughterBranchAddresses(TTree *Tree, std::vector<double> &Momenta, std::vector<double> &MomentaKalmanFit) const;
    /**
     * Function that takes in a sample of double tagged events and determines the raw number of events in each bin and region in beam constrained mass plane
     * @param Tree TTree with double tagged events
     */
    void CalculateBinnedRawYields(TTree *Tree);
    /**
     * Function that calculates the sideband background subtracted yield in bin \f$i\f$
     * @param i Bin number
     */
    double GetBinYield(int i) const;
    /**
     * Function that calculates the total sideband background subtracted yield
     */
    double GetTotalYield() const;
    /**
     * Get the number of events outside phase space
     */
    int GetEventsOutsidePhaseSpace() const;
    /**
     * Get the number of events outside the regions in beam constrained mass space considered
     */
    int GetEventsOutsideMBCSpace() const;
  private:
    /**
     * Number of bins in binning scheme
     */
    int m_NBins;
    /**
     * Phase space parameterisation used to determine which bin events belong to
     */
    AmplitudePhaseSpace m_AmplitudePhaseSpace;
    /**
     * Raw bin yield in the S region, the signal region
     */
    std::vector<int> m_BinYieldS;
    /**
     * Raw bin yield in the A region, real signal event
     */
    std::vector<int> m_BinYieldA;
    /**
     * Raw bin yield in the B region, real tag event
     */
    std::vector<int> m_BinYieldB;
    /**
     * Raw bin yield in the C region, continuum background
     */
    std::vector<int> m_BinYieldC;
    /**
     * Raw bin yield in the D region, combinatorial background
     */
    std::vector<int> m_BinYieldD;
    /**
     * Number of events outside the phase space allowed region
     */
    int m_EventsOutsidePhaseSpace;
    /**
     * Number of events outside the beam constrained mass space considered
     */
    int m_EventsOutsideMBCSpace;
    /**
     * Cuts to apply before counting events
     */
    TCut m_Cuts;
    /**
     * Helper function that calculates the correct bin index
     * Bin numbers are both positive and negative, this function will simply wrap the indexing into a normal array
     */
    int BinIndex(int i) const;
};

#endif
