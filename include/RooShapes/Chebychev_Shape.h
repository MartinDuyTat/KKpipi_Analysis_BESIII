// Martin Duy Tat 18th November 2021
/**
 * This class contains a Chebychev shape of arbitrary order
 * It is used to represent polynomial background in \f$\Delta E\f$ distributions with large combinatorial backgrounds
 */

#ifndef CHEBYCHEV_SHAPE
#define CHEBYCHEV_SHAPE

#include"RooShapes/FitShape.h"

class Chebychev_Shape: public FitShape {
  public:
    /**
     * Default constructor, see FitShape constructor
     */
    Chebychev_Shape(const std::string &Name, const Settings &settings, RooRealVar *x);
    /**
     * Default virtual destructor
     */
    virtual ~Chebychev_Shape() = default;
  private:
    /**
     * Initialize the Chebychev polynomial
     */
    virtual void Initialize();
    /**
     * The order of the polynomial
     */
    int m_Order;
};

#endif
