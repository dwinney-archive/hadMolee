// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef TRIANGLE
#define TRIANGLE

#include "cubature.h"
#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "loop_integrands.hpp"

#include <memory>

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Generic triangle class
    // Contains the masses and call function common to any implementation

    class triangle
    {
        public:

        // Empty, default constructor
        triangle(option mode)
        : _mode(mode)
        {};

        // Parameterized constructor that sets masses
        triangle(option mode, std::array<double,3> external_masses, std::array<double,3> internal_masses)
        : _mode(mode)
        {
            set_external_masses(external_masses);
            set_internal_masses(internal_masses);
        };

        // Evaluate as a function of the initial state invariant mass
        complex eval();
        inline double squared(double s){ return norm(eval()); };

        // Setting functions for the masses
        inline void set_internal_masses(std::array<double,3> m)
        { _imA  = m[0];      _imB  = m[1];      _imC  = m[2]; 
          _imA2 = m[0]*m[0]; _imB2 = m[1]*m[1]; _imC2 = m[2]*m[2]; }
        inline void set_external_masses(std::array<double,3> m)
        { _emA  = m[0];      _emB  = m[1];      _emC  = m[2]; 
          _emA2 = m[0]*m[0]; _emB2 = m[1]*m[1]; _emC2 = m[2]*m[2]; }
        
        // Add a constant width to an internal propagator
        inline void add_width(particle p, double g)
        {
            switch (p) 
            {
                case a: _wA = g; break;
                case b: _wB = g; break;
                case c: _wC = g; break;
                default: return;
            };  
        };

        // Option to change the numerical epsilon used
        void set_ieps(double e){ _eps = e; };

        // Set the maximum number of function calls the integrator is allowed
        void set_max_calls(int n){ _N = n; };

        private:

        // The evaluation depends on what mode is selected and filters through these
        complex nonrelativistic_eval();
        complex relativistic_eval();
        
        // Integrand object which will be used to evaluate the integral with the cubature library
        // Wrapper for the integrand, callable function of feynman parameters
        triangle_integrand integrand = triangle_integrand();
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);

        // Evaluation mode (i.e. relativistic / nonrelativsitic)
        option _mode;

        // External masses
        double _emA,  _emB,  _emC;
        double _emA2, _emB2, _emC2;

        // Internal masses
        double _imA,  _imB,  _imC;
        double _imA2, _imB2, _imC2;
        double _wA = 0., _wB = 0., _wC = 0.; // Widths

        // Default iepsilon perscription
        double _eps = 1.E-4;

        // Number of funciton calls 
        int _N = 1E7;
    };
};

#endif