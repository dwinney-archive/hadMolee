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
    { _imA  = m[0];      _imB  = m[1];      _imC  = m[2]; 
      _imA2 = m[0]*m[0]; _imB2 = m[1]*m[1]; _imC2 = m[2]*m[2]; }
    inline void set_external_masses(array<double,3> m)
    { _emA  = m[0];      _emB  = m[1];      _emC  = m[2]; 
      _emA2 = m[0]*m[0]; _emB2 = m[1]*m[1]; _emC2 = m[2]*m[2]; }

    protected:

    // External masses
    double _emA,  _emB,  _emC;
    double _emA2, _emB2, _emC2;

    // Internal masses
    double _imA,  _imB,  _imC;
    double _imA2, _imB2, _imC2;
};

// ---------------------------------------------------------------------------
// Analytic expression for scalar triangle by treating all particles as nonrelativistic

class nonrelativistic_triangle : public triangle
{
    public:

    // Simple one function and thats it 
    complex<double> eval();
};


// ---------------------------------------------------------------------------
// Relativistic triangle evaluated via Feynman parameters

class relativistic_triangle : public triangle
{
    // -----------------------------------------------------------------------

    public:

    // Evaluate by integrating over Feynman parameters
    complex<double> eval();

    // Option to change the numerical epsilon used
    void set_ieps(double e){ integrand.set_ieps(e); };

    // Set the maximum number of function calls the integrator is allowed
    void set_max_calls(int n){ _N = n; };

    // -----------------------------------------------------------------------

    private: 

    // Number of funciton calls 
    int _N = 1E7;

    // Integrand object which will be used to evaluate the integral with the cubature library
    triangle_integrand integrand;

    // Wrapper for the integrand, callable function of feynman parameters
    static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
};



#endif