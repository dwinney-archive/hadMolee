// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef CONSTANTS
#define CONSTANTS

#include <cmath>
#include <memory>
#include <array>
#include <vector>
#include <iostream>
#include <functional>

#include "complex.hpp"
#include "print.hpp"
#include "debug.hpp"

#ifndef PI
    const double PI       = M_PI;
#endif

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Mathematical constants

    // Fundamental constants
    const double EULER    = 0.57721566490153286060651209008240243104215933593992;
    const double ALPHA    = 1. / 137.;
    const double E        = sqrt(4. * PI * ALPHA);
    
    // Unit conversion factors
    const double DEG2RAD  = (M_PI / 180.);
    const double RAD2DEG  = (180. / M_PI);
    const double NB2GEV   = 2.56819E-6;  // nb -> GeV^{-2} conversion
    const double GEV2NB   = 0.3894E6;    // GeV^{-2} -> nb 

    // Small offsets
    const double  EPS  = 1.E-8;
    const complex IEPS = I*EPS;

    // Three cartesian indexes for ease;
    enum  cartesian_index{x = 0, y = 1, z = 2};
    const cartesian_index C_INDICES[] = { x, y, z };

    // Delta function, made explicit for cartesian_index so you dont need to
    // put cartesian_index::x all the time
    inline bool delta(cartesian_index i, cartesian_index j){ return (i == j); };

    // ---------------------------------------------------------------------------
    // 2022 PDG masses in GeV

    // Light unflavored
    const double M_PION      = 0.13957039;
    const double M_PION0     = 0.1349768;
    const double M_RHO       = 0.77526;

    // Charged charmed
    // const double M_D         = 1.86966;
    // const double M_DSTAR     = 2.01026;
    // const double M_D1        = 2.4221;

    // Use Qians numbers
    const double M_D     = 1.86483;
    const double M_DSTAR = 2.00685;
    const double M_D1    = 2.4208;
    
    // Neutral charmed
    const double M_D0        = 1.86484;
    const double M_DSTAR0    = 2.00685;
    const double M_D10       = 2.412;

    // Charmonia
    const double M_JPSI      = 3.096900;
    const double M_PSI2S     = 3.68610;
    const double M_HC1       = 3.52538;
    const double M_CHIC1     = 3.51067;

    // Exotics
    const double M_Y4260     = 4.2227;
    const double M_ZC3900    = 3.8871;

    // ---------------------------------------------------------------------------
    // 2022 PDG widths in GeV

    const double W_D1       = 31.7E-3; // Charged (narrow) D1
    const double W_D10      = 314.E-3; // Neutral (wide) D1
    const double W_Y4260    = 49.E-3;  
    const double W_ZC3900   = 28.4E-3;

    // ---------------------------------------------------------------------------
    // Couplings associated with HQSS
    
    const double F_PION     = sqrt(2.)*91E-3;  // Pion decay constant in GeV
    const double F_JPSI     = 0.416;

    // D1 -> D* pi
    const double HD_PRIME   = 1.;
    const double H1_D       = sqrt(2./3.)*HD_PRIME/F_PION;
    const double HS_PRIME   = 0.572;
    const double H1_S       = HS_PRIME/sqrt(3.)/F_PION;

    // D* -> D pi
    const double G_TILDE    = 0.57;
    const double G1_PION    = 2.*G_TILDE*sqrt(M_DSTAR*M_D)/F_PION;

    // D* -> D* pi
    const double G2_PION    = 2.*G_TILDE*M_DSTAR/F_PION;

    // Jpsi to charm
    const double G_PSI      = sqrt(M_JPSI) / (2*M_D*F_JPSI);

    const double CT_TRI     = 0.;

    // ---------------------------------------------------------------------------
    // Function for easier comparison of doubles using the EPS value defined above
    // be careful when using this in general purposes since its a fixed-tolerance comparision and not always appropriate

    inline bool is_equal(double a, double b)
    {
        return ( std::abs(a - b) < EPS );
    }

    // Same thing for comparing complex doubles
    inline bool is_equal(complex a, complex b)
    {
        return (is_equal(real(a), real(b)) && is_equal(imag(a), imag(b)));
    };

    // Aliases for special cases of the above
    inline bool is_zero(double a)
    {
        return (std::abs(a) < EPS);
    };

    // ---------------------------------------------------------------------------
    // Kallen triangle function

    // Only way to get a double or int Kallen is if all inputs are double/int
    template <typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };
    // If any of them are complex, return complex
    inline complex Kallen(complex z, double a, double b) { return Kallen<complex>(z, R*a, R*b); };
    inline complex Kallen(double a, complex z, double b) { return Kallen<complex>(R*a, z, R*b); };
    inline complex Kallen(double a, double b, complex z) { return Kallen<complex>(R*a, R*b, z); };
    
    // ---------------------------------------------------------------------------
    // Related to spinors and Lorentz tensor algebra

    // Mostly minus metric
    const double METRIC[4] = {1., -1., -1., -1.};
};

#endif
