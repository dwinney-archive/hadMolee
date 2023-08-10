// Class to evaluate the scalar box function
// We assume external particles are: A + B -> C + D
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "box.hpp"

namespace hadMolee
{
    complex box::eval()
    {
        switch (_mode)
        {
            case (kRelativistic):    return relativistic_eval();
            case (kLoopTools)   :    return looptools_eval();
            default: return std::nan("");
        };
    };

    complex box::looptools_eval()
    {
        // Grab all the masses
        double p01 = _integrand._p01;
        double p02 = _integrand._p02;
        double p03 = _integrand._p03;
        double p12 = _integrand._p12;
        double p13 = _integrand._p13;
        double p23 = _integrand._p23;

        double m0  = _integrand._m0;
        double m1  = _integrand._m1;
        double m2  = _integrand._m2;
        double m3  = _integrand._m3;

        // And widths 
        double w0  = _integrand._w0;
        double w1  = _integrand._w1;
        double w2  = _integrand._w2;
        double w3  = _integrand._w3;
        bool use_complex_masses = is_zero(w0 + w1 + w2 + w3);

        ltini();
        complex integral;
        if ( use_complex_masses ) {integral = D0 (p01, p03, p23, p12, p13, p02, m1, m0, m3, m2);}
        else                      {integral = D0C(p01, p03, p23, p12, p13, p02, m1 - I*csqrt(m1)*w1, m0 - I*csqrt(m0)*w0, m3 - I*csqrt(m3)*w3, m2 - I*csqrt(m2)*w2);}
        ltexi();
        
        return integral / pow(4.*PI, 2.);
    };
    
    complex box::relativistic_eval()
    {
        // Desination for the result and assosiated errors
        double val[3], err[3];

        // Integrate both x and y from 0 to 1
        double min[3] = {0., 0., 0.};
        double max[3] = {1., 1., 1.};


        // TODO: Set relative errors and max calls to actual good values
        // Integrate over x and y
        hcubature(2, wrapped_integrand, &_integrand, 3, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

        // Assemble the result as a complex double
        complex result(val[0], val[1]);
        result /= pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

        return result;
    };

    int box::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
    {
        integrand* Integrand = (integrand *) fdata;

        // Feynman parameters
        double u = in[0], v = in[1], w = in[2];

        double x = u *     v  * (1.-w);
        double y = u * (1.-v) * (1.-w);
        double z = u *              w ;
        double r = 1. - x - y - z; // redundant parameter

        complex result = u*u*(1.-w) * Integrand->eval(r, x, y, z);

        // Split up the real andi imaginary parts to get them out
        fval[0] = std::real(result);
        fval[1] = std::imag(result);

        return 0.;
    };
};