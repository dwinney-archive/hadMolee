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
    integrand.update_masses({_emA2, _emB2, _emC2}, {_imA2, _imB2, _imC2});

    // TODO: Set relative errors and max calls to actual good values
    // Integrate over x and y
    hcubature(2, wrapped_integrand, &integrand, 2, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

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
    double mu_AB = _imA * _imB / (_imA + _imB);
    double mu_BC = _imB * _imC / (_imB + _imC);
    
    // ieps perscription here is to make sure the cuts line up correctly
    complex<double> q_B = sqrt((_emA2 + IEPS)-pow(_emB+_emC, 2.))*sqrt((_emA2 - IEPS)-pow(_emB-_emC, 2.))/(2.*_emA);
    complex<double> E_B =  (_emA2 + _emB2 - _emC2)  / (2.*_emA);

    complex<double> a = pow(q_B*mu_AB/_imA, 2.);

    complex<double> b_AB = _imA + _imB + E_B - _emA;
    complex<double> b_BC = _imB + _imC       - _emA;

    complex<double> c_1 = 2.*mu_BC*b_BC;
    complex<double> c_2 = 2.*mu_AB*b_AB + q_B*q_B*mu_AB/_imA;

    double prefactors = (mu_AB*mu_BC) / (16.*PI*_imA*_imB*_imC);

    complex<double> cut_1 = atan( (c_2 - c_1)        / sqrt( 4.*a*(c_1     - IEPS) ));
    complex<double> cut_2 = atan( (c_2 - c_1 - 2.*a) / sqrt( 4.*a*(c_2 - a - IEPS) ));

    return prefactors / sqrt(a) * (cut_1 - cut_2);
};