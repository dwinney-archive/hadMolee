// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "triangle.hpp"

std::complex<double> hadMolee::relativistic_triangle::eval()
{
    // Desination for the result and assosiated errors
    double val[2], err[2];

    // Integrate both x and y from 0 to 1
    double min[2] = {0., 0.};
    double max[2] = {1., 1.};

    // Fix the "masses" s and t
    integrand.update_masses({_emA2, _emB2, _emC2}, {_imA2, _imB2, _imC2}, {_wA, _wB, _wC});

    // TODO: Set relative errors and max calls to actual good values
    // Integrate over x and y
    hcubature(2, wrapped_integrand, &integrand, 2, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

    // Assemble the result as a complex double
    std::complex<double> result(val[0], val[1]);
    result /= pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

    return result;
};

int hadMolee::relativistic_triangle::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
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

std::complex<double> hadMolee::nonrelativistic_triangle::eval()
{
    std::complex<double> imA = _imA - XI*_wA/2.;
    std::complex<double> imB = _imB - XI*_wB/2.;
    std::complex<double> imC = _imC - XI*_wC/2.;

    std::complex<double> mu_AC = imA * imC / (imA + imC);
    std::complex<double> mu_BC = imB * imC / (imB + imC);
    
    // ieps perscription here is to make sure the cuts line up correctly
    std::complex<double> q_C = sqrt((_emA2 + IEPS)-pow(_emB+_emC, 2.))*sqrt((_emA2 - IEPS)-pow(_emB-_emC, 2.))/(2.*_emA);
    std::complex<double> E_C = sqrt(q_C*q_C + _emC2);

    std::complex<double> a = pow(q_C*mu_AC/_imA, 2.);

    std::complex<double> b_AC = imA + imC + E_C - _emA;
    std::complex<double> b_BC = imB + imC       - _emA;

    std::complex<double> c_1 = 2.*mu_BC*b_BC;
    std::complex<double> c_2 = 2.*mu_AC*b_AC + q_C*q_C*mu_AC/imA;

    std::complex<double> prefactors = (mu_AC*mu_BC) / (16.*PI*imA*imB*imC);

    std::complex<double> cut_1 = atan( (c_2 - c_1)        / sqrt( 4.*a*(c_1     - XI * _eps) ));
    std::complex<double> cut_2 = atan( (c_2 - c_1 - 2.*a) / sqrt( 4.*a*(c_2 - a - XI * _eps) ));

    return prefactors / sqrt(a) * (cut_1 - cut_2);
};