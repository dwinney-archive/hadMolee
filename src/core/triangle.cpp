// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "triangle.hpp"

complex<double> relativistic_triangle::eval()
{
    // Desination for the result and assosiated errors
    double val[2], err[2];

    // Integrate both x and y from 0 to 1
    double min[2] = {0., 0.};
    double max[2] = {1., 1.};

    // Fix the "masses" s and t
    integrand.update_masses({_ema2, _emb2, _emc2}, {_ima2, _imb2, _imc2});

    // TODO: Set relative errors and max calls to actual good values
    // Integrate over x and y
    hcubature(2, wrapped_integrand, &integrand, 2, min, max, 2E7, 0, 1e-2, ERROR_INDIVIDUAL, val, err);

    // Assemble the result as a complex double
    std::complex<double> result(val[0], val[1]);
    result *= 2.;                  // Factor of 2 from the normalization of dF_3 integration measure
    result /= 2. * pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

    return result;
};

int relativistic_triangle::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
  triangle_integrand* integrand = (triangle_integrand *) fdata;

  // Feynman parameters
  double x = in[0] * in[1];
  double y = in[0] * (1. - in[1]);
  double z = 1. - x - y;

  std::complex<double> result = in[0] * integrand->eval(x, y, z);

  // Split up the real andi imaginary parts to get them out
  fval[0] = std::real(result);
  fval[1] = std::imag(result);

  return 0.;
};

complex<double> nonrelativistic_triangle::eval()
{
    complex<double> c1 = 2.*mu12()*b12();
    complex<double> c2 = 2.*mu23()*b12() + 2.*qB()*qB()*mu23()/_imc;
    complex<double> a = pow(qB()*mu23()/_imc, 2.);
    complex<double> N = mu12()*mu23() / (16.*PI*_ima*_imb*_imc);

    complex<double> term1 = atan( XR*(c2 - c1)       / sqrt(4.*a*(c1 - IEPS))     );
    complex<double> term2 = atan( XR*(c2 - c1- 2.*a) / sqrt(4.*a*(c2 - a - IEPS)) );

    // debug("qb", qB());
    // debug("mu23", mu23());
    // debug("c1", c1);
    // debug("c2", c2);
    // debug("a", a);
    // debug("term1", term1);
    // debug("term2", term2);

    return N / sqrt(XR*a)*(term1 - term2);
};