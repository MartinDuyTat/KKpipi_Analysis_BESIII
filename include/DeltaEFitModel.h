// Martin Duy Tat 17th November 2021
/**
 * This class sets up the fit model for \f$\Delta E\f$
 */

#ifndef DELTAEFITMODEL
#define DELTAEFITMODEL

#include<vector>
#include"RooArgList.h"
#include"RooAddPdf.h"
#include"RooRealVar.h"
#include"Settings.h"
#include"RooShapes/FitShape.h"

class DeltaEFitModel {
  public:
    /**
     * Constructor that sets up the fit model
     * @param settings The settings object containing all the model parameters
     * @param x The independent fit variable
     */
    DeltaEFitModel(const Settings &settings, RooRealVar *x);
    /**
     * Destructor that destroys all PDF components
     */
    ~DeltaEFitModel();
    /**
     * Get the full \f$\Delta E\f$ fit model
     */
    RooAbsPdf* GetModel() const;
    /**
     * Get the PDF component
     */
    RooAbsPdf* GetModelComponent(int i) const;
  private:
    /**
     * Vector of all PDF components
     */
    std::vector<FitShape*> m_ModelComponents;
    /**
     * List of all PDFs
     */
    RooArgList m_ModelPDFs;
    /**
     * The full fit model of \f$\Delta E\f$
     */
    RooAddPdf* m_FullModel = nullptr;
    /**
     * Fit settings
     */
    Settings m_Settings;
    /**
     * The independent fit variable
     */
    RooRealVar *m_x;
    /**
     * Initialize the signal component
     */
    void InitializeSignal();
    /**
     * Initialize the combinatorial background component
     */
    void InitializeCombinatorial();
    /**
     * Initialize full fit shape
     */
    void InitializeFullShape();
};

#endif
