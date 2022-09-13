// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef TRIANGLE
#define TRIANGLE

#include "cubature.h"
#include "constants.hpp"
#include "loop_integrands.hpp"

// ---------------------------------------------------------------------------
// Generic triangle class
// Contains the masses and call function common to any implementation

class triangle
{
    public:

    triangle()
    {};

    triangle(array<double,3> external_masses, array<double,3> internal_masses)
    {
        set_external_masses(external_masses);
        set_internal_masses(internal_masses);
    };

    // Evaluate as a function of the initial state invariant mass
    virtual complex<double> eval() = 0;
    double squared(double s){ return norm(eval()); };

    // Setting functions for the masses
    inline void set_internal_masses(array<double,3> m)
    { _ima  = m[0];      _imb  = m[1];      _imc  = m[2]; 
      _ima2 = m[0]*m[0]; _imb2 = m[1]*m[1]; _imc2 = m[2]*m[2]; }
    inline void set_external_masses(array<double,3> m)
    { _ema  = m[0];      _emb  = m[1];      _emc  = m[2]; 
      _ema2 = m[0]*m[0]; _emb2 = m[1]*m[1]; _emc2 = m[2]*m[2]; }

    protected:

    // External masses
    double _ema,  _emb,  _emc;
    double _ema2, _emb2, _emc2;

    // Internal masses
    double _ima,  _imb,  _imc;
    double _ima2, _imb2, _imc2;
};

// ---------------------------------------------------------------------------
// Analytic expression for scalar triangle by treating all particles as nonrelativistic

class nonrelativistic_triangle : public triangle
{
    public:

    // Empty constructor for setting masses later
    nonrelativistic_triangle()
    : triangle()
    {};

    // Else specify the internal masses
    nonrelativistic_triangle(array<double,3> internal_masses)
    : triangle({0., 0., 0.}, internal_masses)
    {};

    complex<double> eval();

    private:

    // Reduced masses
    inline complex<double> mu12(){ return _ima * _imb / (_ima + _imb); };
    inline complex<double> mu23(){ return _imb * _imc / (_imb + _imc); };

    // Mass differences
    inline complex<double> b12() { return _ima + _imb - _ema; };
    inline complex<double> b23() { return _imb + _imc + EB() - _ema; };

    // Energy and momentum of external particle b in the rest frame of particle a
    inline complex<double> qB() { return sqrt(Kallen(XR*_ema2, XR*_emb2, XR*_emc2)) / (2. * _ema); };
    inline complex<double> EB() { return sqrt(qB()*qB() + _emb2); };
};


// ---------------------------------------------------------------------------
// Relativistic triangle evaluated via Feynman parameters

class relativistic_triangle : public triangle
{
    // -----------------------------------------------------------------------

    public:

    // Empty constructor to set masses later
    relativistic_triangle()
    : triangle()
    {};

    // Parameterized constructor with all masses
    relativistic_triangle(array<double,3> external_masses, array<double,3> internal_masses)
    : triangle(external_masses, internal_masses)
    {}

    // Evaluate by integrating over Feynman parameters
    complex<double> eval();

    // -----------------------------------------------------------------------

    private: 

    // Integrand object which will be used to evaluate the integral with the cubature library
    triangle_integrand integrand;

    // Wrapper for the integrand, callable function of feynman parameters
    static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
};



#endif