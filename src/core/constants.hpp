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
#include <complex>
#include <memory>
#include <array>
#include <vector>
#include <iostream>
#include "debug.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Mathematical constants

    // Generic particle labels
    enum particle{a, b, c, d};

    // We use complex numbers a lot so we define this shortened data type
    using complex = std::complex<double>;

    // Fundamental constants
    const double PI       = M_PI;
    const double EULER    = 0.57721566490153286060651209008240243104215933593992;
    const double ALPHA    = 1. / 137.;
    const double E        = sqrt(4. * PI * ALPHA);

    // Unit complex and real values
    const complex XR  (1., 0.);
    const complex XI  (0., 1.);
    
    // Unit conversion factors
    const double DEG2RAD  = (M_PI / 180.);
    const double RAD2DEG  = (180. / M_PI);
    const double NB2GEV   = 2.56819E-6;  // nb -> GeV^{-2} conversion
    const double GEV2NB   = 0.3894E6;    // GeV^{-2} -> nb 

    // Small offsets
    const double                EPS  = 1.E-8;
    const complex IEPS  = XI*EPS;

    // Three cartesian indexes for ease;
    enum  cartesian_index{x = 0, y = 1, z = 2};
    const cartesian_index C_INDICES[] = { x, y, z };

    // Delta function, made explicit for cartesian_index so you dont need to
    // put cartesian_index::x all the time
    inline bool delta(cartesian_index i, cartesian_index j){ return (i == j); };

    // Many objects have different evaluations which can be controlled through this option flag
    enum option{ relativistic, nonrelativistic, LoopTools };

    // ---------------------------------------------------------------------------
    // 2022 PDG masses in GeV

    // Light unflavored
    const double M_PION      = 0.13957039;
    const double M_PION0     = 0.1349768;

    // Charged charmed
    const double M_D         = 1.86966;
    const double M_DSTAR     = 2.01026;
    const double M_D1        = 2.4221;

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

    const double W_D1       = 31.3E-3; // Charged (narrow) D1
    const double W_D10      = 314.E-3; // Neutral (wide) D1
    const double W_Y4260    = 49.E-3;  
    const double W_ZC3900   = 28.4E-3;

    // ---------------------------------------------------------------------------
    // Couplings associated with HQSS and previous analyses 

    // Qin & Qiang (2016) [arXiv:1509.01398]
    const double MY_QQ2016        = 4.217;    // Pole mass
    const double YNM_WIDTH_QQ2016 = 0.056;    // Non-molecular width of the Y-meson
    const double F_Y_QQ2016       = 1./0.063; // Y-meson decay constant (unitless)
    const double YBARE_QQ2016     = 10.88;    // YD1D coupling in GeV-1/2  
    const double ZBARE_QQ2016     = 0.77;     // ZD*D coupling in GeV-1/2
    const double G1_QQ2016        = 29.50;    // Elastic D1D  coupling in GeV-2
    const double A_QQ2016         = 12.67;    // Overall S-wave normalization coeff in GeV-5/2
    const double B_QQ2016         = -15.23;   // Const S-wave polynomail coeff in GeV2

    // Cleven et al (2014) 	[arXiv:1310.2190]
    const double A_CLEVEN2014     = 6.72;     // Overall S-wave normalization coeff in GeV -1
    const double B_CLEVEN2014     = -15.28;   // Const S-wave polynomial coeff in GeV2 

    // HQSS constants
    // Primed couplings have 
    const double HPRIME_UPPER     = 1.31;     // Upperbound of D1 D* pi coupling using full width (Gamma = 25 MeV) in GeV-1
    const double F_PION           = 132.E-3;  // Pion decay constant in GeV

    // ---------------------------------------------------------------------------
    // Overload of multiplcation for a bool and complex<double>

    inline complex operator * (const bool & a, const complex & b)
    {
        if (a)
            return b;
        else
            return complex(0, 0);
    };
    inline complex operator * (const complex & b, const bool & a)
    {
        if (a)
            return b;
        else
            return complex(0, 0);
    };

    // Same thing but with int
    inline complex operator * (const int & a, const complex & b)
    {
        if (a != 0) return complex(double(a) * real(b), double(a) * imag(b));
        else        return 0.*XR;
    };
    inline complex operator * (const complex & b, const int & a)
    {
        if (a != 0) return complex(double(a) * real(b), double(a) * imag(b));
        else        return 0.*XR;
    };

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
        return (abs(a) < EPS);
    };

    // ---------------------------------------------------------------------------
    // Helper function to help "static cast" derived unique_ptr to base unique_ptr
    // fron this https://stackoverflow.com/questions/36120424/alternatives-of-static-pointer-cast-for-unique-ptr

    template<typename TO, typename FROM>
    std::unique_ptr<TO> static_unique_pointer_cast(std::unique_ptr<FROM>&& old)
    {
        return std::unique_ptr<TO>{static_cast<TO*>(old.release())};
    };

    // ---------------------------------------------------------------------------
    // Kallen triangle function
    template <typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    // ---------------------------------------------------------------------------
    // Related to spinors and Lorentz tensor algebra

    // Mostly minus metric
    const double METRIC[4] = {1., -1., -1., -1.};

    // Gamma matrix vector in Dirac basis
    const complex GAMMA[4][4][4] =
    {
        //gamma0
        { { 1., 0., 0., 0. },
        { 0., 1., 0., 0. },
        { 0., 0., -1., 0. },
        { 0., 0., 0., -1. } },
        //gamma1
        { { 0., 0., 0., 1. },
        { 0., 0., 1., 0. },
        { 0., -1., 0., 0. },
        { -1., 0., 0., 0. } },
        //gamma2
        { { 0., 0., 0., -XI },
        { 0., 0., XI, 0. },
        { 0.,  XI, 0., 0. },
        { -XI, 0., 0., 0. } },
        //gamma3
        { { 0., 0., 1., 0. },
        { 0., 0., 0., -1. },
        { -1., 0., 0., 0. },
        { 0., 1., 0., 0. } }
    };

    // Gamma_5
    const complex GAMMA_5[4][4] =
    {
        { 0., 0., 1., 0. },
        { 0., 0., 0., 1. },
        { 1., 0., 0., 0. },
        { 0., 1., 0., 0. }
    };
};

#endif
