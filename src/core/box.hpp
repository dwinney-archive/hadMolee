// Class to evaluate the scalar box function
// We assume external particles are: A + B -> C + D
// and use the notation for particle masses as in [1]
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/1005.2076
// ---------------------------------------------------------------------------


#ifndef BOX_HPP
#define BOX_HPP

#include "constants.hpp"
#include "cubature.h"
#include "clooptools.h"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Generic box class
    // Contains the masses and call function common to any implementation

    class box
    {
        public:

        // Empty, default constructor
        box(int mode)
        : _mode(mode)
        {
            if (mode != kRelativistic)
            {
                warning("box", "Initialized with unavailable option, numerical output may be NaN...");
            }
        };

        // Parameterized constructor that sets masses
        box(int mode, std::array<double,4> external_masses, std::array<double,4> internal_masses)
        : _mode(mode)
        {
            set_external_masses(external_masses);
            set_internal_masses(internal_masses);
        };

        // Evaluate as a function of the initial state invariant mass
        complex eval();
        inline complex eval(std::array<double,4> external, std::array<double,4> internal, double s, double t)
        {
            set_external_masses(external); 
            set_internal_masses(internal);
            set_invariant_masses(s, t);
            return eval();
        };
        inline double squared(double s){ return norm(eval()); };

        // Setting function for the masses of the internal masses in the loop
        // Order of arguments starts at top left leg and goes counter-clockwise
        inline void set_external_masses(std::array<double,4> m)
        {  _integrand._p01 = m[0]*m[0]; _integrand._p12 = m[1]*m[1]; _integrand._p23 = m[2]*m[2]; _integrand._p03 = m[3]*m[3]; }

        // Setting function for the (on-shell) masses of external particles
        // Order of arguments starts at top, horizontal propagator and goes counter-clockwise
        inline void set_internal_masses(std::array<double,4> m)
        { _integrand._m0  = m[0]*m[0]; _integrand._m1  = m[1]*m[1]; _integrand._m2  = m[2]*m[2]; _integrand._m3  = m[3]*m[3]; }

        // Invariant masses of the 02 and 13 systems 
        // Equivalently particles AC and AB respectively
        inline void set_invariant_masses(double s, double t)
        {  
            _integrand._p02    = s;  _integrand._p13    = t; 
        };
        
        // Add a constant width to an internal propagator
        inline void add_width(int i, double w)
        {
            switch (i)
            {
                case 0: _integrand._w0 = w; return;
                case 1: _integrand._w1 = w; return;
                case 2: _integrand._w2 = w; return;
                case 3: _integrand._w3 = w; return;
                default: return;
            };
        };

        // Option to change the numerical epsilon used
        void set_ieps(double e){ _integrand._eps = e; };

        // Set the maximum number of function calls the integrator is allowed
        void set_max_calls(int n){ _N = n; };

        // Different options for evaluation
        static const int kRelativistic    = 0;
        static const int kLoopTools       = 1;

        // ---------------------------------------------------------------------------

        private:

        // Flag whether or not we're using brute-forced box or LoopTools
        int _mode;
    
        // ---------------------------------------------------------------------------
        // These are related to the box evaluated with cubature integration over Feynman parameters

        // Evaluate the box by integrating over Feynman parameters
        complex relativistic_eval();

        // Container class to help us sneak non-static member data into hcubature
        struct integrand
        {
            // ---------------------------------------------------------------------------
            // Saves mass / energy variables 

            // These are all masses SQUARED!! (GeV^2)!!!!
            double _p01, _p02, _p03, _p12, _p13, _p23; 
            double _m0, _m1, _m2, _m3;
            double _w0, _w1, _w2, _w3;

            // Default epsilon
            double _eps = EPS;

            // Evaluate the integrand at fixed values of the feynman parameters
            inline complex eval(double x0, double x1, double x2, double x3)
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
        };

        // Integrand object which will be used to evaluate the integral with the cubature library
        integrand _integrand;

        // Wrapper for the integrand, callable function of feynman parameters
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
        
        // Maximum number of function calls allowed
        int _N = 1E7;

        // ----------------------------------------------------------------------
        // Methods for evaluating using LoopTools
        complex looptools_eval();
    };
};

#endif