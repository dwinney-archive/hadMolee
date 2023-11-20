// Specific implementations of amplitudes relevant for the D* D pi final state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DSDPI_D1_HPP
#define DSDPI_D1_HPP

#include <memory>

#include "kinematics.hpp"
#include "breit_wigner.hpp"
#include "amplitude.hpp"
#include "constants.hpp"
#include "triangle.hpp"
#include "molecule.hpp"
#include "Y(4260).hpp"

namespace hadMolee::DsDpi
{
    class D1 : public amplitude_base
    {
        public: 

        D1(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "generic_D1")
        : amplitude_base(key, xkinem, V, 0, "generic_D1", id),
          _Y(get_molecular_component(_V))
          {};

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _A * D1_vertex(i, j);
        };

        // Options or evaluating only one piece
        static const int kSwaveOnly = 1;
        static const int kDwaveOnly = 2;

        protected:

        inline complex D1_vertex(cartesian_index i, cartesian_index j)
        {
            double gs = (_option == kDwaveOnly) ? 0. : H1_S;
            double gd = (_option == kSwaveOnly) ? 0. : H1_D;

            double swave = gs*sqrt(_mpc*_mpc + _mc2)            * delta(i,j);
            double dwave = gd*_mpc*_mpc*(3.*phat_1(i)*phat_1(j) - delta(i,j));

            return sqrt(M_D1*M_DSTAR)*(swave + dwave);
        };

        // Energy dependent (but scalar) coupling
        complex _A = 0;
        
        // Need direct access to Y molecular nature
        molecule _Y;  
    };

    // ---------------------------------------------------------------------------
    // Transition involving D1DsD triangle with the Z as well as 
    // No free parameters since this assume to be known background

    class triangle : public D1
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        triangle(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_triangle")
        : D1(key, xkinem, V, id),
          _T(hadMolee::triangle::kNonrelativistic), 
          _Zc(make_molecule(M_D, M_DSTAR))
        {
            // Set up Zc propagator
            _Zc->set_parameters({3.9, 50.E-3, 4.66/sqrt(M_DSTAR*M_D*3.9)});

            // Set up triangle
            _T.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _T.add_width(2, W_D1); 
        };

        // -----------------------------------------------------------------------
        private:

        inline void recalculate()
        {
            // Update arguments with floating Y and Z meson masses for the triangle
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _A = -2.0*_T.eval();

            // Gather all couplings
            _A *= _Y->coupling()/sqrt(2.) * sqrt(M_D1*M_D*_Y->mass());  
            // Skip D-wave D1 -> D* pi which we factor out
            _A *= _Zc->coupling()         * sqrt(M_D*M_DSTAR*_Zc->mass()); 
            _A *= _Zc->propagator(_sab);
            _A *= _Zc->coupling()         * sqrt(M_D*M_DSTAR*_Zc->mass()); 
        };

        // Scalar triangle function
        hadMolee::triangle _T;

        // This channel recieves contribution from the Z(3900) & Y (molecular) couplings
        molecule _Zc;
    };

    // ---------------------------------------------------------------------------
    // Also include the D1 tree transition 
    // No free parameters since this assume to be known background

    class tree : public D1
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        tree(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_tree")
        : D1(key, xkinem, V, id),
          _D1(breit_wigner::kNonrelativistic, M_D1, W_D1)
        {};

        // -----------------------------------------------------------------------
        private:

        inline void recalculate()
        {            
            _A  = _Y->coupling() / sqrt(2.) * sqrt(_Y->mass()*M_D1*M_D);  
            _A *= _D1.eval(_sac); // D1 propagator
            // Skip D1 -> D* pi which we factor out and put above
        };

        // In addition we have the tree level transition
        breit_wigner _D1;
    };
};

#endif