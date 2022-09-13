// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "triangle.hpp"

complex<double> relativistic_triangle::eval()
{
    auto F = [&](double x)
    {
        complex<double> iEPS = _eps*XI;
        complex<double> a, b, c, d, yp, ym;
        a = _ema2 + iEPS;
        b = - (_imb2 - _imc2) - (_ema2 + iEPS)*(x-1.) + (_emb2 - _emc2)*x;
        c = _ima2*x + _imc2*(1.-x) - _emb2*(1.-x)*x;

        d = b*b - 4.*a*c;

        yp = (- b + sqrt(XR*d))/ (2.*a);
        ym = (- b - sqrt(XR*d))/ (2.*a);

        return (log(1. - x + ym) - log(1. - x + yp) - log(ym) + log(yp)) / sqrt(XR*d);
    };

    complex<double> result = gauss_kronrod<double, 61>::integrate(F, 0., 1., 0., 1.E-6, NULL);
    return result;
};

complex<double> nonrelativistic_triangle::eval()
{
    complex<double> c1 = 2.*mu12()*b12();
    complex<double> c2 = 2.*mu23()*b12() + 2.*qB()*qB()*mu23()/_imc;
    complex<double> a = pow(qB()*mu23()/_imc, 2.);
    complex<double> N = mu12()*mu23() / (16.*PI*_ima*_imb*_imc);

    complex<double> term1 = atan( XR*(c2 - c1)       / sqrt(4.*a*(c1 - IEPS))     );
    complex<double> term2 = atan( XR*(c2 - c1- 2.*a) / sqrt(4.*a*(c2 - a - IEPS)) );

    debug("qb", qB());
    debug("mu23", mu23());
    debug("c1", c1);
    debug("c2", c2);
    debug("a", a);
    debug("term1", term1);
    debug("term2", term2);

    return N / sqrt(XR*a)*(term1 - term2);
};