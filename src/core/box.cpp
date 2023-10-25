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
            case (kRelativistic):   return relativistic_eval();
            case (kLoopTools)   :   return looptools_eval();
            case (kTOPT)        :   return TOPT_eval();
            default: return std::nan("");
        };
    };

    // ---------------------------------------------------------------------------
    // LoopTools stuff

    std::array<complex,4> box::eval_vector()
    {
        if (_mode != kLoopTools) return {std::nan(""), std::nan(""), std::nan(""), std::nan("")};

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
        double w0  = _integrand._w0 + 1E-6;
        double w1  = _integrand._w1 + 1E-6;
        double w2  = _integrand._w2 + 1E-6;
        double w3  = _integrand._w3 + 1E-6;

        ltini();
        setwarndigits(60);
        std::array<int,4> dd = {dd0, dd1, dd2, dd3};
        std::array<complex,4> vector;
        for (int i = 0; i < 4; i++)
        {
           vector[i] = D0iC(dd[i], p01, p03, p23, p12, p13, p02, m1 - I*csqrt(m1)*w1, m0 - I*csqrt(m0)*w0, m3 - I*csqrt(m3)*w3, m2 - I*csqrt(m2)*w2) / pow(4.*PI, 2.);
        };
        ltexi();
        
        return vector;
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
        setwarndigits(60);
        complex integral;
        if ( use_complex_masses ) {integral = D0i (dd0, p01, p03, p23, p12, p13, p02, m1, m0, m3, m2);}
        else                      {integral = D0iC(dd0, p01, p03, p23, p12, p13, p02, m1 - I*csqrt(m1)*w1, m0 - I*csqrt(m0)*w0, m3 - I*csqrt(m3)*w3, m2 - I*csqrt(m2)*w2);}
        ltexi();
        
        return integral / pow(4.*PI, 2.);
    };
    
    // ---------------------------------------------------------------------------
    // Brute force relativistic stuff

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

        return 0;
    };

    complex box::integrand::eval(double x0, double x1, double x2, double x3)
    {
        complex Y01, Y02, Y03, Y12, Y13, Y23;

        complex M0 = _m0 - I*csqrt(_m0)*_w0 + IEPS;
        complex M1 = _m1 - I*csqrt(_m1)*_w1 + IEPS;
        complex M2 = _m2 - I*csqrt(_m2)*_w2 + IEPS;
        complex M3 = _m3 - I*csqrt(_m3)*_w3 + IEPS;

        Y01 = M0 + M1 - _p01;
        Y02 = M0 + M2 - _p02;
        Y03 = M0 + M3 - _p03;
        Y12 = M1 + M2 - _p12;
        Y13 = M1 + M3 - _p13;
        Y23 = M2 + M3 - _p23;        

        complex D = x0*x0*M0 + x1*x1*M1   + x2*x2*M2  + x3*x3*M3
                                + x0*x1*Y01  + x0*x2*Y02 + x0*x3*Y03
                                            + x1*x2*Y12 + x1*x3*Y13 
                                                        + x2*x3*Y23;
        return pow(D - I*_eps, -2.);
    };

    // ---------------------------------------------------------------------------
    // TOPT stuff

    // We integrate over the loop 3-momentum
    complex box::TOPT_eval()
    {
        // Desination for the result and assosiated errors
        double val[3], err[3];

        // Integrate over {r, theta, phi} 
        // The r integration mapped to r = tan(pi x/2) where x goes from 0 to 1
        double min[3] = {0.,   0,  0.   };
        double max[3] = {1.,  PI,  2.*PI};

        // Integrate over 3D momentum
        hcubature(2, wrapped_integrand_TOPT, &_integrand, 3, min, max, _N, 0, 1e-6, ERROR_INDIVIDUAL, val, err);

        // Assemble the result as a complex double
        complex result(val[0], val[1]);
        result /= pow(2.*PI, 3.); // Divide by the phase-space normalization 

        return result;
    };

    int box::wrapped_integrand_TOPT(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
    {
        integrand* Integrand = (integrand *) fdata;

        // r is the integration momentum
        double x     = in[0];
        double theta = in[1];
        double phi   = in[2];
        double r     = tan(PI*x/2.);

        // Jacobian includes the r -> x change as well as the angular element
        double jacobian = (PI/2.) * r*r / pow(cos(PI*x/2.), 2.) * sin(theta);

        // Relegate all evaluation to the integrand
        complex result = jacobian * Integrand->eval_TOPT(r, theta, phi);

        // Split up the real andi imaginary parts to get them out
        fval[0] = std::real(result);
        fval[1] = std::imag(result);

        return 0;
    };

    // At this stage we asssume that all external variables have been set
    complex box::integrand::eval_TOPT(double r, double theta, double phi)
    {
        // The total energy is the decay particle mass
        complex E = sqrt(_p01) + I*_eps;

        // Dot products of 3 momenta
        complex r_dot_p03   = r*_mom03* cos(theta);
        complex r_dot_p23   = r*_mom23*(cos(theta)*_cos02 + cos(phi)*sin(theta)*sqrt(1.-_cos02*_cos02));
        complex p03_dot_p23 = _mom03*_mom23*_cos02; 

        // Calculate omega which are the energies of the internal particles
        std::array<complex,6> omega;
        omega[0] = csqrt(_m0 + I*(_w0 + _eps)/2. + r*r);
        omega[1] = csqrt(_m1 + I*(_w1 + _eps)/2. + r*r);
        omega[2] = csqrt(_m2 + I*(_w2 + _eps)/2. + r*r + _mom03*_mom03 + _mom23*_mom23 + 2.*(p03_dot_p23 - r_dot_p03 - r_dot_p23));
        omega[3] = csqrt(_m3 + I*(_w3 + _eps)/2. + r*r + _mom03*_mom03 - 2.*r_dot_p03);

        std::array<complex,3> G;

        // First time ordening
        G[0] = E - omega[0] - omega[1];
        G[1] = E - _E03     - omega[3] - omega[1];
        G[2] = E - _E03     - omega[3] - omega[2] - _E12;

        if (is_zero(abs(G[0]))) print("G0 is zero!");
        if (is_zero(abs(G[1]))) print("G1 is zero!");
        if (is_zero(abs(G[2]))) print("G2 is zero!");

        return 1./(16.*omega[0]*omega[1]*omega[2]*omega[3]*G[0]*G[1]*G[2]);


        // // Propagators first index is time ordering second index is cut
        // std::array< std::array<complex,3>, 2> G;

        // // Second time ordering
        // G[1][0] = E - omega[0] - omega[1];
        // G[1][1] = E - omega[0] - omega[2] - _E12;
        // G[1][2] = E - _E03     - omega[3] - omega[2] - _E12;
        
        // complex sum_of_orderings = 0.;
        // for (auto time_ordering : G)
        // {
        //     complex product_of_Gs = 1.;
        //     for (auto propagator : time_ordering)
        //     {
        //         product_of_Gs *= propagator;
        //     }
        //     sum_of_orderings += 1./product_of_Gs;
        // };

        // return sum_of_orderings / (16.*omega[0]*omega[1]*omega[2]*omega[3]);
    };
};