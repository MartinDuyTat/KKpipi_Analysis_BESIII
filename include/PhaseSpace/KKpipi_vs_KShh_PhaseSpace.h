// Martin Duy Tat 20th December 2021
/**
 * KKpipi_vs_KShh_PhaseSpace inherits from KKpipi_PhaseSpace, but saves the two Dalitz coordinates of the KShh phase space
 * Only absolute bin number on the tag side is determined
 */

#ifndef KKPIPI_VS_KSHH_PHASESPACE
#define KKPIPI_VS_KSHH_PHASESPACE

#include<string>
#include"TTree.h"
#include"TH2F.h"
#include"TLorentzVector.h"
#include"KKpipi_PhaseSpace.h"

class KKpipi_vs_KShh_PhaseSpace: public KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that saves the tag mode, Dalitz coordinates and binning scheme
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     * @param ReconstructedBins Set to true to calculate the reconstructed bins
     * @param TrueBins Set to true to calculate the true bins
     * @param Mode "KSpipi" or "KSKK"
     */
    KKpipi_vs_KShh_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, const std::string &Mode);
    /**
     * Get the correct phase space bin with flavour tag
     */
    virtual std::pair<int, int> Bin() const;
    /**
     * Get the correct true phase space bin with flavour tag
     */
    virtual std::pair<int, int> TrueBin();
  private:
    /**
     * Read bin number from lookup table
     */
    int LookUpBinNumber(double M2Plus, double M2Minus) const;
    /**
     * Map a Dalitz point outside of phase space back inside the boundary
     */
    int GetMappedKShhBin(double M2Plus, double M2Minus) const;
    /**
     * Get bin number
     * If Dalitz point is outside phase space it's automatically mapped back inside the boundary
     */
    int GetKShhBin() const;
    /**
     * Get true bin number
     * If Dalitz point is outside phase space it's automatically mapped back inside the boundary
     */
    int GetTrueKShhBin();
    /**
     * The KShh mode
     */
    const std::string m_Mode;
    /**
     * Flag that is 1 if Kalman fit worked
     */
    int m_KalmanFitSuccess;
    /**
     * The four-momentum of \f$K_S\f$
     */
    TLorentzVector m_KS_P;
    /**
     * The four-momentum of \f$h^+\f$
     */
    TLorentzVector m_hPlus_P;
    /**
     * The four-momentum of \f$h^-\f$
     */
    TLorentzVector m_hMinus_P;
    /**
     * The four-momentum of \f$K_S\f$ after Kalman fit
     */
    TLorentzVector m_KS_P_KalmanFit;
    /**
     * The four-momentum of \f$h^+\f$ after Kalman fit
     */
    TLorentzVector m_hPlus_P_KalmanFit;
    /**
     * The four-momentum of \f$h^-\f$ after Kalman fit
     */
    TLorentzVector m_hMinus_P_KalmanFit;
    /**
     * The KShh binning scheme
     */
    TH2F *m_BinningScheme;
    /**
     * Set KSpipi branch addresses
     */
    void SetKSpipiBranchAddresses(TTree *Tree);
    /**
     * Set branch addresses common to both KSpipi and KSKK
     */
    void SetKShhBranchAddresses(TTree *Tree);
};

#endif
