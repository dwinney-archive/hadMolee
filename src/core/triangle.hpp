// Class to evaluate the scalar triangle function
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "cubature.h"
#include "constants.hpp"
#include "clooptools.h"

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
        triangle(int mode)
        : _mode(mode)
        {};

        // Parameterized constructor that sets masses
        triangle(int mode, std::array<double,3> external_masses, std::array<double,3> internal_masses)
        : _mode(mode)
        {
            set_external_masses(external_masses);
            set_internal_masses(internal_masses);
        };

        // Evaluate as a function of the initial state invariant mass
        complex eval();

        // Evaluate the vector decomposition
        std::array<complex,3> eval_vector();

        // Return the function squared
        inline double squared(double s){ return norm(eval()); };

        // Setting functions for the masses
        inline void set_internal_masses(std::array<double,3> m)
        { 
            _integrand._ima  = m[0];      _integrand._imb  = m[1];      _integrand._imc  = m[2]; 
            _integrand._ima2 = m[0]*m[0]; _integrand._imb2 = m[1]*m[1]; _integrand._imc2 = m[2]*m[2];

            // If internal masses change, reset all the widths
            _integrand._wa = 0.; _integrand._wb = 0.; _integrand._wc = 0.; 
        }
        inline void set_external_masses(std::array<double,3> m)
        {
            _integrand._ema  = m[0];      _integrand._emb  = m[1];      _integrand._emc  = m[2]; 
            _integrand._ema2 = m[0]*m[0]; _integrand._emb2 = m[1]*m[1]; _integrand._emc2 = m[2]*m[2];
        }
        
        // Add a constant width to an internal propagator
        inline void add_width(int i, double g)
        {
            switch (i) 
            {
                case 1: _integrand._wa = g; break;
                case 2: _integrand._wb = g; break;
                case 3: _integrand._wc = g; break;
                default: return;
            };  
        };

        // Option to change the numerical epsilon used
        void set_ieps(double e){ _integrand._eps = e; };

        // Set the maximum number of function calls the integrator is allowed
        void set_max_calls(int n){ _N = n; };

        // Different options for evaluation
        static const int kRelativistic    = 0;
        static const int kNonrelativistic = 1;
        static const int kLoopTools       = 2;

        private:
          // Auxillary class to help interfacing with the cubature library

        struct integrand
        {
            // Evaluate the integrand at fixed values of the feynman parameters
            inline std::complex<double> eval(double x1, double x2, double x3)
            {
                complex ima2 = _ima2 - I*sqrt(_ima2)*_wa;
                complex imb2 = _imb2 - I*sqrt(_imb2)*_wb;
                complex imc2 = _imc2 - I*sqrt(_imc2)*_wc;
                
                std::complex<double> D =  x1*ima2
                                        + x2*imb2
                                        + x3*imc2
                                        - x1*x3*_emb2
                                        - x1*x2*_emc2 
                                        - x2*x3*_ema2;
                return 1. / (D - I*_eps);
            };

            // Default epsilon
            double _eps = EPS;

            // Need to be able to access the masses
            double _ema,  _emb,  _emc;
            double _ima,  _imb,  _imc;

            double _ema2, _emb2, _emc2;
            double _ima2, _imb2, _imc2;
            double _wa = EPS, _wb = EPS, _wc = EPS;
        };

        // The evaluation depends on what mode is selected and filters through these
        complex nonrelativistic_eval();
        complex relativistic_eval();
        complex looptools_eval();
        
        // Integrand object which will be used to evaluate the integral with the cubature library
        // Wrapper for the integrand, callable function of feynman parameters
        integrand _integrand;
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);

        // Evaluation mode (i.e. relativistic / nonrelativsitic)
        int _mode;

        // Number of funciton calls 
        int _N = 1E7;
    };
};

#endif