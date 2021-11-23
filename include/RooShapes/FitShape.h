// Martin Duy Tat 17th November 2021
/**
 * This is an abstract class that contains a general fit shape component
 * Any class that inherits from this one needs to initialize all the parameters and the PDF by overloading the Initialize() function
 */

#ifndef FITSHAPE
#define FITSHAPE

#include<map>
#include<string>
#include"RooAbsPdf.h"
#include"RooRealVar.h"
#include"Settings.h"

class FitShape {
  public:
    /**
     * Constructor that saves the settings
     * @param Name Name of this PDF
     * @param settings Settings describing this PDF
     * @param x The independent variable in the fit
     */
    FitShape(const std::string &Name, const Settings &settings, RooRealVar *x);
    /**
     * Empty virtual destructor
     */
    virtual ~FitShape();
    /**
     * Get a pointer to the PDF
     */
    RooAbsPdf* GetPDF();
    /**
     * Get a pointer to the yield variable
     */
    RooAbsReal* GetYield();
    /**
     * Use this function to change the yield parameterisation from an independent yield variable to a fixed ratio relative to another yield variable
     * @param SignalYield Signal yield variable that we want to fix the yield to
     * @param BackgroundToSignalYieldRatio The product of this and the signal yield gives us the yield of this component
     */
    void UseRelativeYield(RooRealVar *SignalYield, double BackgroundToSignalYieldRatio);
  protected:
    /**
     * Settings for the PDF shape
     */
    Settings m_Settings;
    /**
     * Name of this component
     */
    std::string m_Name;
    /**
     * Map of RooRealVar objects that define the PDF, labelled by their names
     */
    std::map<std::string, RooRealVar*> m_Parameters;
    /**
     * The actual PDF
     */
    RooAbsPdf *m_PDF = nullptr;
    /**
     * Pointer to the independent variable in this fit
     */
    RooRealVar *m_x;
  private:
    /**
     * Initialize all the parameters and the PDF
     */
    virtual void Initialize() = 0;
    /**
     * RooRealVar containing the yield of this component
     */
    RooAbsReal *m_Yield = nullptr;
};

#endif
