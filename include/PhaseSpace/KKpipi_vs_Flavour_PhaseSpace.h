// Martin Duy Tat 26th November 2021
/**
 * KKpipi_vs_Flavour_PhaseSpace inherits from KKpipi_PhaseSpace, but also saves the charge of the kaon on the tag side
 * The kaon charge gives us information about the flavour, so we can swap the sign of the bin number if necessary
 */

#ifndef KKPIPI_VS_FLAVOUR_PHASESPACE
#define KKPIPI_VS_FLAVOUR_PHASESPACE

#include<utility>
#include"PhaseSpace/KKpipi_PhaseSpace.h"
#include"TTree.h"

class KKpipi_vs_Flavour_PhaseSpace: public KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that saves the branch of the charge of the tag kaon
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     */
    KKpipi_vs_Flavour_PhaseSpace(TTree *Tree, int Bins);
    /**
     * Get the correct phase space bin with flavour tag
     */
    virtual std::pair<int, int> Bin() const;
  private:
    /**
     * Charge of the tag kaon
     */
    double m_KaonCharge;
};

#endif
