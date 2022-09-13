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

// Dressed propagator for a hadronic molecule
// This assumes a meson X is an S-wave molecule of constituent hadrons m1 and m2
class hadronic_propagator : public propagator
{
    public:

    // Constructor requires three masses:
    // the (dressed) pole mass of the resonance and the masse of the two constituent mesons 
    // Optional bool changes the behavior of the self-energy correction
    // If true use the simpler MS-bar renormalized construction
    // if false default to using the full 2-point loop function
    hadronic_propagator(double pole_mass, double m1, double m2, bool use_msbar = true)
    : propagator(pole_mass), _m1(m1), _m2(m2), _msbar(use_msbar)
    {};

    // Set the three free parameters
    inline void set_params(std::array<double, 3> par)
    {
        _binding      = par[0];
        _elastic      = par[1];
        _nonmol_width = par[2];
    };

    // Assemble the self-energy and width into the dressed non-relativistic propagator
    complex<double> eval(double E);

    // Function specifying the self-energy correction brough upon from the m1 m2 channel
    complex<double> self_energy(double E); 

    // Allow one of the particles, m1, to have a constant width
    inline void add_constitent_width(double gam){ _width1 = gam; };

    protected:

    // If true use the simpler MS-bar renormalized construction
    // if false default to using the full 2-point loop function
    double _msbar = true;

    // Constituent meson masses
    double _m1, _m2;
    double _width1 = 0.;
    inline double reduced_mass(){ return _m1*_m2 / (_m1 + _m2); };

    // Three free parameters
    double _binding;  // Coupling constant of X -> m1 m2 interaction
    double _elastic;  // Coupling for elastic m1 m2 -> m1 m2 scattering
    double _nonmol_width;             // Constant other decays than into m1 m2
};

#endif