// Abstract class for an amplitude of a particles scalar propagator
// as well as some implementations for the Y and Z mesons
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef PROPAGATORS
#define PROPAGATORS

#include "constants.hpp"

#include <boost/math/differentiation/finite_difference.hpp>

// First the abstract class that can be called from anywhere
class propagator
{
    public: 

    // Default constructor only requires a 'bare' mass 
    propagator(double mass)
    : _M(mass)
    {};

    // Main function, evaluates as a function of one variable
    virtual complex<double> eval(double x) = 0;
    inline double squared(double x){ return norm(eval(x)); };

    protected:

    // Pole mass is the only 
    double _M;
};

// Relativistic BW with constant width
class relativistic_BW : public propagator
{
    public: 
    // Default constructor takes in a pole mass and a constant width
    relativistic_BW(double mass, double width)
    : propagator(mass), _Gamma(width)
    {};

    inline complex<double> eval(double s)
    {
        complex<double> D =  (s - _M*_M) + XI * _M*_Gamma;
        return XI / D;
    };

    private:

    // Constant width
    double _Gamma;
};

// Non-relativistic BW with constant width
class nonrelativistic_BW : public propagator
{
    public: 
    // Default constructor takes in a pole mass and a constant width
    nonrelativistic_BW(double mass, double width)
    : propagator(mass), _Gamma(width)
    {};

    inline complex<double> eval(double E)
    {
        complex<double> D =  (E - _M) + XI * _Gamma/2.;
        return XI / (2. * D);
    };

    private:

    // Constant width
    double _Gamma;
};

#endif