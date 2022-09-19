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
const double M_Y4260     = 4.2227;
const double M_Z3900     = 3.8871;

// ---------------------------------------------------------------------------
// 2022 PDG widths in GeV

const double W_D1       = 31.3E-3; // Charged (narrow) D1
const double W_D10      = 314.E-3; // Neutral (wide) D1
const double W_Y4260    = 49.E-3;  
const double W_Z3900    = 28.4E-3;

// ---------------------------------------------------------------------------
// Couplings associated with HQSS and previous analyses 

const double C_Z        = 0.77;     // GeV^-2 [from arXiv:1310.2190 ]
const double F_PION     = 92.1E-3;  // Pion decay constant [from arXiv:2201.08253 ]

// The couplings below get the D-wave coupling from the D2 decay then require rest of D1 width to come from S-wave
// [from arXiv:2001.05884]
const double HP_D        = 1.17;     // GeV^-2 D-wave D1 Ds pi coupling
const double HP_S        = 0.57;     //        S-wave D1 Ds pi coupling

// Here we assume the full D1 width is saturated by the D-wave 
const double HP          = 2.21394;     // GeV^-2 D-wave D1Dspi [from arXiv:2001.05884]


// ---------------------------------------------------------------------------
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
const std::complex<double> GAMMA[4][4][4] =
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
const complex<double> GAMMA_5[4][4] =
{
    { 0., 0., 1., 0. },
    { 0., 0., 0., 1. },
    { 1., 0., 0., 0. },
    { 0., 1., 0., 0. }
};

// ---------------------------------------------------------------------------

#endif
