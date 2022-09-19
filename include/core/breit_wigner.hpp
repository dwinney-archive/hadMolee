// Abstract class for an amplitude of a particles scalar propagator
// as well as some implementations for the Y and Z mesons
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef BREITWIGNER
#define BREITWIGNER

#include "constants.hpp"

// First the abstract class that can be called from anywhere
class breit_wigner
{
    public: 

    // Default constructor only requires a 'bare' mass 
    breit_wigner(double mass, double width = 0.)
    : _mass(mass), _width(width)
    {};

    // Main function, evaluates as a function of one variable
    virtual complex<double> eval(double x) = 0;
    inline double squared(double x){ return norm(eval(x)); };

    protected:

    // Pole mass is only really needed parameter 
    double _mass;

    // But may also allow a width
    double _width;
};

// Relativistic BW with constant width
class relativistic_BW : public breit_wigner
{
    public: 
    // Default constructor takes in a pole mass and a constant width
    relativistic_BW(double mass, double width)
    : breit_wigner(mass, width)
    {};

    inline complex<double> eval(double s)
    {
        complex<double> D =  (s - _mass*_mass) + XI * _mass*_width + IEPS;
        return XI / D;
    };
};

// Non-relativistic BW with constant width
class nonrelativistic_BW : public breit_wigner
{
    public: 
    // Default constructor takes in a pole mass and a constant width
    nonrelativistic_BW(double mass, double width)
    : breit_wigner(mass, width)
    {};

    inline complex<double> eval(double E)
    {
        complex<double> D =  (E - _mass) + XI * _width/2. + IEPS;
        return XI / (2. * D);
    };
};

#endif