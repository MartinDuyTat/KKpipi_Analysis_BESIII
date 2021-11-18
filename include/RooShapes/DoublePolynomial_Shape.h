// Martin Duy Tat 17th November 2021
/**
 * This class contains a piece-wise continuous second-order polynomial, centered around zero
 * It is used to represent polynomial background in \f$\Delta E\f$ distributions with large combinatorial backgrounds
 */

#ifndef DOUBLEPOLYNOMIAL_SHAPE
#define DOUBLEPOLYNOMIAL_SHAPE

#include"RooShapes/FitShape.h"

class DoublePolynomial_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    DoublePolynomial_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
  private:
    /**
     * Initialize a double polynomial
     */
    virtual void Initialize();
};

#endif
