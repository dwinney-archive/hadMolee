// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "triangle.hpp"

namespace hadMolee
{
    complex triangle::eval()
    {
        switch (_mode)
        {
            case (kRelativistic):    return relativistic_eval();
            case (kNonrelativistic): return nonrelativistic_eval();
            case (kLoopTools):       return looptools_eval();
            default: return std::nan("");
        };
    };
    
    complex triangle::looptools_eval()
    {
        // Grab all the masses
        double p1 = _integrand._ema2;
        double p2 = _integrand._emb2;
        double p3 = _integrand._emc2;

        double m1  = _integrand._ima2;
        double m2  = _integrand._imb2;
        double m3  = _integrand._imc2;

        // And widths 
        double w1  = _integrand._wa + 1.E-6;
        double w2  = _integrand._wb + 1.E-6;
        double w3  = _integrand._wc + 1.E-6;

        ltini();
        setwarndigits(60);
        complex integral = C0C(p1, p2, p3, m2 - I*csqrt(m2)*w2, m3 - I*csqrt(m3)*w3, m1 - I*csqrt(m1)*w1);
        ltexi();
        
        return - integral / pow(4.*PI, 2.);
    };

    std::array<complex,3> triangle::eval_vector()
    {
        if (_mode != kLoopTools) return {nan(""), nan(""), nan("")};

        // Grab all the masses
        double p1 = _integrand._ema2;
        double p2 = _integrand._emb2;
        double p3 = _integrand._emc2;

        double m1  = _integrand._ima2;
        double m2  = _integrand._imb2;
        double m3  = _integrand._imc2;

        // And widths 
        double w1  = _integrand._wa + 1.E-6;
        double w2  = _integrand._wb + 1.E-6;
        double w3  = _integrand._wc + 1.E-6;

        ltini();
        setwarndigits(60);
        std::array<int,3> cc = {cc0, cc1, cc2};
        std::array<complex,3> vector;
        for (int i = 0; i < 3; i++)
        {
            vector[i] = - C0iC(cc[i], p1, p2, p3, m2 - I*csqrt(m2)*w2, m3 - I*csqrt(m3)*w3, m1 - I*csqrt(m1)*w1) / pow(4.*PI, 2.);

        }
        ltexi();
        
        return vector;
    };

    complex triangle::relativistic_eval()
    {
        // Desination for the result and assosiated errors
        double val[2], err[2];

        // Integrate both x and y from 0 to 1
        double min[2] = {0., 0.};
        double max[2] = {1., 1.};

        // TODO: Set relative errors and max calls to actual good values
        // Integrate over x and y
        hcubature(2, wrapped_integrand, &_integrand, 2, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

        // Assemble the result as a complex double
        complex result(val[0], val[1]);
        result /= pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

        return result;
    };

    int triangle::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
    {
        integrand* Integrand = (integrand *) fdata;

        // Feynman parameters
        double x = in[0] * in[1];
        double y = in[0] * (1. - in[1]);
        double z = 1. - x - y;

        complex result = in[0] * Integrand->eval(x, y, z);

        // Split up the real andi imaginary parts to get them out
        fval[0] = std::real(result);
        fval[1] = std::imag(result);

        return 0;
    };

    complex triangle::nonrelativistic_eval()
    {
        complex imA = _integrand._ima - I*_integrand._wa/2.;
        complex imB = _integrand._imb - I*_integrand._wb/2.;
        complex imC = _integrand._imc - I*_integrand._wc/2.;

        complex mu_AC = imA * imC / (imA + imC);
        complex mu_BC = imB * imC / (imB + imC);
        
        // ieps perscription here is to make sure the cuts line up correctly
        complex q_C = sqrt((_integrand._ema2 + IEPS)-pow(_integrand._emb+_integrand._emc, 2.))*sqrt((_integrand._ema2 - IEPS)-pow(_integrand._emb-_integrand._emc, 2.))/(2.*_integrand._ema);
        complex E_C = sqrt(q_C*q_C + _integrand._emc2);

        complex a = pow(q_C*mu_AC/_integrand._ima, 2.);

        complex b_AC = imA + imC + E_C - _integrand._ema;
        complex b_BC = imB + imC       - _integrand._ema;

        complex c_1 = 2.*mu_BC*b_BC;
        complex c_2 = 2.*mu_AC*b_AC + q_C*q_C*mu_AC/imA;

        complex prefactors = (mu_AC*mu_BC) / (16.*PI*imA*imB*imC);

        complex cut_1 = atan( (c_2 - c_1)        / sqrt( 4.*a*(c_1     - I * _integrand._eps) ));
        complex cut_2 = atan( (c_2 - c_1 - 2.*a) / sqrt( 4.*a*(c_2 - a - I * _integrand._eps) ));
        
        return prefactors / sqrt(a) * (cut_1 - cut_2);
    };
};