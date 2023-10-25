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
            if (mode != kRelativistic && mode != kLoopTools && mode != kTOPT)
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

        // Output the components of the vector integral
        std::array<complex,4> eval_vector();
        inline std::array<complex,4> eval_vector(std::array<double,4> external, std::array<double,4> internal, double s, double t)
        {
            set_external_masses(external); 
            set_internal_masses(internal);
            set_invariant_masses(s, t);
            return eval_vector();
        };

        // Setting function for the masses of the internal masses in the loop
        // Order of arguments starts at top left leg and goes counter-clockwise
        inline void set_external_masses(std::array<double,4> m)
        {  
            _integrand._p01 = m[0]*m[0]; _integrand._p12 = m[1]*m[1]; _integrand._p23 = m[2]*m[2]; _integrand._p03 = m[3]*m[3]; 
            if (_mode == kTOPT) _integrand.update_TOPT();
        }

        // Setting function for the (on-shell) masses of external particles
        // Order of arguments starts at top, horizontal propagator and goes counter-clockwise
        inline void set_internal_masses(std::array<double,4> m)
        { 
            _integrand._m0  = m[0]*m[0]; _integrand._m1  = m[1]*m[1]; _integrand._m2  = m[2]*m[2]; _integrand._m3  = m[3]*m[3];
            // if internal masses change clear all the widths
            _integrand._w0 = 0.; _integrand._w1 = 0.; _integrand._w2 = 0.; _integrand._w3 = 0.;   
        }

        // Invariant masses of the 02 and 13 systems 
        // Equivalently particles AC and AB respectively
        inline void set_invariant_masses(double s, double t)
        {  
            _integrand._p02    = s;  _integrand._p13    = t; 
            if (_mode == kTOPT) _integrand.update_TOPT();
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
        static const int kTOPT            = 2;

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
            complex eval(double x0, double x1, double x2, double x3);

            // Evaluate using the TOPT expression
            complex eval_TOPT(double r, double theta, double phi);

            // 3 momenta of each external particle
            inline void update_TOPT()
            {
                // Assume partilce 01 rest-frame defines "lab frame"
                // In the particle labeling of hadMolee::particle this is V
                // we also thus have _m0 = s is the total invariant energy

                // Each other particle's three momentum
                // Remmeber _m and _p's are already squared
                _mom12 = csqrt(Kallen(_m0, _p02, _p12)) / (2.*sqrt(_m0));
                _mom03 = csqrt(Kallen(_m0, _p13, _p03)) / (2.*sqrt(_m0));
                _mom23 = csqrt(Kallen(_m0, _p01 + _p12 + _p23 + _p03 - _p01 - _p13, _p23)) / (2.*sqrt(_m0));

                // Energies of externals
                _E12 = csqrt( _mom12*_mom12 + _p12);
                _E03 = csqrt( _mom03*_mom03 + _p03);
                _E23 = csqrt( _mom23*_mom23 + _p23);

                // Relative angle between 03 and 23
                _cos02 = (2.*_E03*_E23 - _p02 + _p03 + _p23) / (2.*_mom03*_mom23);

                line();
                print("Recalculating...");
                print("mom12", _mom12);
                print("mom03", _mom03);
                print("mom23", _mom23);

                print("E12", _E12);
                print("E03", _E03);
                print("E23", _E23);

                print("cos", _cos02);
            };

            complex  _E12, _E23, _E03;       // Energies 
            complex  _mom12, _mom23, _mom03; // Momenta
            complex  _cos02;                  // Relative angle between 03 and 23
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

        // ----------------------------------------------------------------------
        // Methods for evaluating using TOPT expression
        complex TOPT_eval();

        static int wrapped_integrand_TOPT(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
    };
};

#endif