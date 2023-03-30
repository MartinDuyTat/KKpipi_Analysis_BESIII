// Martin Duy Tat 17th November 2021
/**
 * This class contains a double Gaussian RooFit shape where the mean and width are parameterised as a ratio
 */

#ifndef DOUBLEGAUSSIANRATIO_SHAPE
#define DOUBLEGAUSSIANRATIO_SHAPE

#include<string>
#include"RooShapes/FitShape.h"
#include"Settings.h"

class DoubleGaussianRatio_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    DoubleGaussianRatio_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
    /**
     * Default virtual destructor
     */
    virtual ~DoubleGaussianRatio_Shape() = default;
  private:
    /**
     * Initialize a double Gaussian
     */
    virtual void Initialize();
};

#endif
