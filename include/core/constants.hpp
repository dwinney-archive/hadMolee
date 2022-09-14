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
#include "debug.hpp"

// ---------------------------------------------------------------------------
// Mathematical constants

const double PI       = M_PI;
const double DEG2RAD  = (M_PI / 180.);
const double RAD2EG   = (180. / M_PI);
const double ALPHA    = 1. / 137.;
const double E        = sqrt(4. * PI * ALPHA);
const double EULER    = 0.57721566490153286060651209008240243104215933593992;

// Unit complex and real values
const complex<double> XR  (1., 0.);
const complex<double> XI  (0., 1.);

// Small offsets
const double           EPS  = 1.E-6;
const complex<double> IEPS  = XI*EPS;

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
const double M_Y         = 4.2227;
const double M_Z         = 3.8871;

// Overload of multiplcation for a bool and complex<double>
inline complex<double> operator * (const bool & a, const complex<double> & b){
    if (a != 0)
        return b;
     else
        return complex<double>(0, 0);
};
inline complex<double> operator * (const complex<double> & b, const bool & a){
    if (a != 0)
        return b;
     else
        return complex<double>(0, 0);
};

// Kallen triangle function
template <typename T>
inline T Kallen(T x, T y, T z)
{
    return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
};


// ---------------------------------------------------------------------------

#endif
