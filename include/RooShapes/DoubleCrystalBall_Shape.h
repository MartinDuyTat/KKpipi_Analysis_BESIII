// Martin Duy Tat 17th November 2021
/**
 * This class contains a double Crystal Ball RooFit shape
 */

#ifndef DOUBLECRYSTALBALL_SHAPE
#define DOUBLECRYSTALBALL_SHAPE

#include<string>
#include"RooShapes/FitShape.h"
#include"Settings.h"

class DoubleCrystalBall_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    DoubleCrystalBall_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
    /**
     * Default virtual destructor
     */
    virtual ~DoubleCrystalBall_Shape() = default;
  private:
    /**
     * Initialize a double Gaussian
     */
    virtual void Initialize();
};

#endif
