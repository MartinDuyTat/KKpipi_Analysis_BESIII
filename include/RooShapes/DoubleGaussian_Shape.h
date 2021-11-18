// Martin Duy Tat 17th November 2021
/**
 * This class contains a double Gaussian RooFit shape
 */

#ifndef DOUBLEGAUSSIAN_SHAPE
#define DOUBLEGAUSSIAN_SHAPE

#include<string>
#include"RooShapes/FitShape.h"
#include"Settings.h"

class DoubleGaussian_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    DoubleGaussian_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
  private:
    /**
     * Initialize a double Gaussian
     */
    virtual void Initialize();
};

#endif
