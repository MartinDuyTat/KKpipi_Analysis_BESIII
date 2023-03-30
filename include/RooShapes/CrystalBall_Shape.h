// Martin Duy Tat 17th November 2021
/**
 * This class contains a Crystal Ball RooFit shape
 */

#ifndef CRYSTALBALL_SHAPE
#define CRYSTALBALL_SHAPE

#include<string>
#include"RooShapes/FitShape.h"
#include"Settings.h"

class CrystalBall_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    CrystalBall_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
    /**
     * Default virtual destructor
     */
    virtual ~CrystalBall_Shape() = default;
  private:
    /**
     * Initialize a double Gaussian
     */
    virtual void Initialize();
};

#endif
