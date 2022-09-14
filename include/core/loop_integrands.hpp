// Auxillary classes used to evaluate the Feynman parameter integrands for loop integrals
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef INTEGRANDS
#define INTEGRANDS

#include "constants.hpp"

// Auxillary class to help interfacing with the cubature library
class triangle_integrand
{
    // -------------------------------------------------------------------

    public:

    // Empty constructor cus we dont need anything
    triangle_integrand(){};

    // Get masses squared from the parent triangle function
    inline void update_masses(array<double,3> ex_m2, array<double,3> in_m2)
    {
        _ema2 = ex_m2[0]; _emb2 = ex_m2[1]; _emc2 = ex_m2[2];
        _ima2 = in_m2[0]; _imb2 = in_m2[1]; _imc2 = in_m2[2];
    };

    // Evaluate the integrand at fixed values of the feynman parameters
    inline complex<double> eval(double x1, double x2, double x3)
    {
        double D = x1*_ima2 + x2*_imb2 + x3*_imc2 - x1*x3*_emb2 - x1*x2*_emc2 - x2*x3*_ema2;
        return 1. / (D - XI*_eps);
    };

    // Set the numerical iepsilon used
    inline void set_ieps(double e){ _eps = e; };

    // -------------------------------------------------------------------

    private:

    // Default epsilon
    double _eps = EPS;

    // Need to be able to access the masses
    double _ema2, _emb2, _emc2;
    double _ima2, _imb2, _imc2;
};

#endif