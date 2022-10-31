// Class to evaluate the scalar box function
// We assume external particles are: A + B -> C + D
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef BOX
#define BOX

#include "cubature.h"
#include "constants.hpp"
#include "loop_integrands.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Generic box class
    // Contains the masses and call function common to any implementation

    class box
    {
        public:

        // Empty, default constructor
        box()
        {};

        // Parameterized constructor that sets masses
        box(std::array<double,4> external_masses, std::array<double,4> internal_masses)
        {
            set_external_masses(external_masses);
            set_internal_masses(internal_masses);
        };

        // Nothing special needed for destructor
        virtual ~box() = default;

        // Evaluate as a function of the initial state invariant mass
        virtual complex eval() = 0;
        double squared(double s){ return norm(eval()); };

        // Setting functions for the masses
        inline void set_internal_masses(std::array<double,4> m)
        { _imA  = m[0];      _imB  = m[1];      _imC  = m[2];      _imD  = m[3];
        _imA2 = m[0]*m[0]; _imB2 = m[1]*m[1]; _imC2 = m[2]*m[2]; _imD2 = m[3]*m[3]; }
        inline void set_external_masses(std::array<double,4> m)
        { _emA  = m[0];      _emB  = m[1];      _emC  = m[2];      _emD  = m[3];
        _emA2 = m[0]*m[0]; _emB2 = m[1]*m[1]; _emC2 = m[2]*m[2]; _emD2 = m[3]*m[3]; }
        inline void set_invariant_masses(double s, double t){ _s = s; _t = t; }
        
        // Add a constant width to an internal propagator
        inline void add_width_a(double g){ _wA = g; };
        inline void add_width_b(double g){ _wB = g; };
        inline void add_width_c(double g){ _wC = g; };
        inline void add_width_d(double g){ _wD = g; };

        protected:

        // External masses
        double _emA,  _emB,  _emC,  _emD;
        double _emA2, _emB2, _emC2, _emD2;

        // Internal masses
        double _imA,  _imB,  _imC,  _imD;
        double _imA2, _imB2, _imC2, _imD2;
        double _wA = 0., _wB = 0., _wC = 0., _wD = 0.; // Widths

        // We also need two invariant masses
        double _s;  // Assumed to be the CD system
        double _t;  // Assumed to be the BD system
    };


    // ---------------------------------------------------------------------------
    // Relativistic box evaluated via explicit numerical integration of Feynman parameters

    class relativistic_box : public box
    {
        // -----------------------------------------------------------------------

        public:

        relativistic_box()
        : box()
        {};

        relativistic_box(std::array<double,4> external_masses, std::array<double,4> internal_masses)
        : box(external_masses, internal_masses)
        {};


        // Evaluate by integrating over Feynman parameters
        complex eval();

        // Option to change the numerical epsilon used
        void set_ieps(double e){ integrand.set_ieps(e); };

        // Set the maximum number of function calls the integrator is allowed
        void set_max_calls(int n){ _N = n; };

        // -----------------------------------------------------------------------

        private: 

        // Number of funciton calls 
        int _N = 1E7;

        // Integrand object which will be used to evaluate the integral with the cubature library
        box_integrand integrand;

        // Wrapper for the integrand, callable function of feynman parameters
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
    };
};

#endif