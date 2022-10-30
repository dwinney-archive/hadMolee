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
#include <memory>

namespace hadMolee
{
    // Currentl availale options for the triangle implementation
    enum triangle_options{relativistic, nonrelativistic, LoopTools};
    class triangle;
    class nonrelativistic_triangle;
    class relativistic_triangle;

    // ---------------------------------------------------------------------------
    // Generic triangle class
    // Contains the masses and call function common to any implementation

    class triangle
    {
        public:

        // Empty, default constructor
        triangle()
        {};

        // Parameterized constructor that sets masses
        triangle(std::array<double,3> external_masses, std::array<double,3> internal_masses)
        {
            set_external_masses(external_masses);
            set_internal_masses(internal_masses);
        };

        // Nothing special needed for destructor
        virtual ~triangle() = default;

        // Evaluate as a function of the initial state invariant mass
        virtual std::complex<double> eval() = 0;
        double squared(double s){ return norm(eval()); };

        // Setting functions for the masses
        inline void set_internal_masses(std::array<double,3> m)
        { _imA  = m[0];      _imB  = m[1];      _imC  = m[2]; 
          _imA2 = m[0]*m[0]; _imB2 = m[1]*m[1]; _imC2 = m[2]*m[2]; }
        inline void set_external_masses(std::array<double,3> m)
        { _emA  = m[0];      _emB  = m[1];      _emC  = m[2]; 
          _emA2 = m[0]*m[0]; _emB2 = m[1]*m[1]; _emC2 = m[2]*m[2]; }
        
        // Add a constant width to an internal propagator
        inline void add_width_a(double g){ _wA = g; };
        inline void add_width_b(double g){ _wB = g; };
        inline void add_width_c(double g){ _wC = g; };

        protected:

        // External masses
        double _emA,  _emB,  _emC;
        double _emA2, _emB2, _emC2;

        // Internal masses
        double _imA,  _imB,  _imC;
        double _imA2, _imB2, _imC2;
        double _wA = 0., _wB = 0., _wC = 0.; // Widths
    };

    // ---------------------------------------------------------------------------
    // Analytic expression for scalar triangle by treating all particles as nonrelativistic

    class nonrelativistic_triangle : public triangle
    {
        public:

        nonrelativistic_triangle()
        : triangle()
        {};

        nonrelativistic_triangle(std::array<double,3> external_masses, std::array<double,3> internal_masses)
        : triangle(external_masses, internal_masses)
        {};

        // Simple 
        std::complex<double> eval();

        // Option to change the numerical epsilon used
        void set_ieps(double e){ _eps = e; };

        private:

        // Default iepsilon perscription
        double _eps = EPS;
    };


    // ---------------------------------------------------------------------------
    // Relativistic triangle evaluated via Feynman parameters

    class relativistic_triangle : public triangle
    {
        // -----------------------------------------------------------------------

        public:

        relativistic_triangle()
        : triangle()
        {};

        relativistic_triangle(std::array<double,3> external_masses, std::array<double,3> internal_masses)
        : triangle(external_masses, internal_masses)
        {};


        // Evaluate by integrating over Feynman parameters
        std::complex<double> eval();

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

    inline std::unique_ptr<triangle> make_triangle(triangle_options x = nonrelativistic)
    {
        switch (x)
        {
            case nonrelativistic: 
            {
                std::unique_ptr<triangle> T = std::make_unique<nonrelativistic_triangle>();
                return T;
            }
            case relativistic: 
            {
                std::unique_ptr<triangle> T = std::make_unique<relativistic_triangle>();
                return T;
                break;
            }
            default: return nullptr;
        }
        return nullptr;
    };
};

#endif