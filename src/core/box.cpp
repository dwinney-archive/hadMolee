// Class to evaluate the scalar box function
// We assume external particles are: A + B -> C + D
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "box.hpp"


std::complex<double> hadMolee::relativistic_box::eval()
{
    // Desination for the result and assosiated errors
    double val[3], err[3];

    // Integrate both x and y from 0 to 1
    double min[3] = {0., 0., 0.};
    double max[3] = {1., 1., 1.};

    // Fix the "masses", s and ,t
    integrand.update_masses({_s, _t}, {_emA2, _emB2, _emC2, _emD2}, {_imA2, _imB2, _imC2, _imD2}, {_wA, _wB, _wC, _wD});

    // TODO: Set relative errors and max calls to actual good values
    // Integrate over x and y
    hcubature(2, wrapped_integrand, &integrand, 3, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

    // Assemble the result as a complex double
    std::complex<double> result(val[0], val[1]);
    result /= pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

    return result;
};

int hadMolee::relativistic_box::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
  box_integrand* integrand = (box_integrand *) fdata;

  // Feynman parameters
  double u = in[0], v = in[1], w = in[2];

  double x = u *     v  * (1.-w);
  double y = u * (1.-v) * (1.-w);
  double z = u *              w ;
  double r = 1. - x - y - z; // redundant parameter

  std::complex<double> result = u*u*(1.-w) * integrand->eval(r, x, y, z);

  // Split up the real andi imaginary parts to get them out
  fval[0] = std::real(result);
  fval[1] = std::imag(result);

  return 0.;
};