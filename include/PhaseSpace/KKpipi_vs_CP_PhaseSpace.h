// Martin Duy Tat 21st December 2021
/**
 * KKpipi_vs_CP_PhaseSpace inherits from KKpipi_PhaseSpace
 * Since the tag side is a CP tag, the signal conjugate bins are equivalent to the regular bins
 */

#ifndef KKPIPI_VS_CP_PHASESPACE
#define KKPIPI_VS_CP_PHASESPACE

#include<utility>
#include"PhaseSpace/KKpipi_PhaseSpace.h"
#include"TTree.h"

class KKpipi_vs_CP_PhaseSpace: public KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that saves the branch of the charge of the tag kaon
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     * @param ReconstructedBins Set to true to calculate the reconstructed bins
     * @param TrueBins Set to true to calculate the true bins
     * @param KSKK_binning Set to true to determine the true \f$K_SKK\f$ bins
     */
    KKpipi_vs_CP_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, bool KSKK_binning = false);
    /**
     * Get the correct phase space bin with CP tag
     */
    virtual std::pair<int, int> Bin() const;
    /**
     * Get the correct true phase space bin with CP tag
     */
    virtual std::pair<int, int> TrueBin();
};

#endif
