// Martin Duy Tat 20th December 2021
/**
 * KKpipi_vs_K0hh_PhaseSpace inherits from KKpipi_PhaseSpace, but saves the two Dalitz coordinates of the K0hh phase space
 * Only absolute bin number on the tag side is determined
 */

#ifndef KKPIPI_VS_K0HH_PHASESPACE
#define KKPIPI_VS_K0HH_PHASESPACE

#include<string>
#include"TTree.h"
#include"TH2F.h"
#include"TLorentzVector.h"
#include"KKpipi_PhaseSpace.h"

class KKpipi_vs_K0hh_PhaseSpace: public KKpipi_PhaseSpace {
  public:
    /**
     * Constructor that saves the tag mode, Dalitz coordinates and binning scheme
     * @param Tree TTree with double tag events
     * @param Bins Number of bins in KKpipi phase space binning
     * @param ReconstructedBins Set to true to calculate the reconstructed bins
     * @param TrueBins Set to true to calculate the true bins
     * @param Mode "KSpipi", "KSKK", "KLpipi", "KLKK"
     * @param KSKK_binning Set to true to determine the true \f$K_SKK\f$ bins
     * @param KKpipiPartReco Set to true if the KKpipi mode is partially reconstructed
     * @param KStoKLBackground Set to true if we're looking at KS background in KL tags
     */
    KKpipi_vs_K0hh_PhaseSpace(TTree *Tree, int Bins, bool ReconstructedBins, bool TrueBins, const std::string &Mode, bool KSKK_binning = false, bool KKpipiPartReco = false, bool KStoKLBackground = false);
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
     * Get bin number
     * If Dalitz point is outside phase space it's automatically mapped back inside the boundary
     */
    int GetK0hhBin() const;
    /**
     * Get true bin number
     * If Dalitz point is outside phase space it's automatically mapped back inside the boundary
     */
    int GetTrueK0hhBin();
    /**
     * The K0 mode, either "KS" or "KL"
     */
    const std::string m_K0Mode;
    /**
     * The hh mode, either "pipi" or "KK"
     */
    const std::string m_hhMode;
    /**
     * Flag that is 1 if Kalman fit worked
     */
    int m_KalmanFitSuccess;
    /**
     * The four-momentum of \f$K_0\f$
     */
    TLorentzVector m_K0_P;
    /**
     * The four-momentum of \f$h^+\f$
     */
    TLorentzVector m_hPlus_P;
    /**
     * The four-momentum of \f$h^-\f$
     */
    TLorentzVector m_hMinus_P;
    /**
     * The four-momentum of \f$K_0\f$ after Kalman fit
     */
    TLorentzVector m_K0_P_KalmanFit;
    /**
     * The four-momentum of \f$h^+\f$ after Kalman fit
     */
    TLorentzVector m_hPlus_P_KalmanFit;
    /**
     * The four-momentum of \f$h^-\f$ after Kalman fit
     */
    TLorentzVector m_hMinus_P_KalmanFit;
    /**
     * The K0hh binning scheme
     */
    TH2F *m_BinningScheme;
    /**
     * Set K0hh branch addresses
     */
    void SetK0hhBranchAddresses(TTree *Tree);
    /**
     * Flag that is true if we're considering KS to KL background
     */
    bool m_KStoKLBackground;
};

#endif
