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
    // ---------------------------------------------------------------------------
    // Transition involving D1DsD triangle with the Z as well as 
    // No free parameters since this assume to be known background

    class triangle : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        triangle(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_triangle")
        : amplitude_base(key, xkinem, V, 0, "DsDpi_dwave", id),
          _T(hadMolee::triangle::kNonrelativistic), 
          _Zc(make_molecule(M_D, M_DSTAR)),
          _Y(get_molecular_component(_V))
        {
            // Set up Zc propagator
            _Zc->set_parameters({3.9, 50.E-3, 4.66/sqrt(M_DSTAR*M_D*3.9)});

            // Set up triangle
            _T.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _T.add_width(2, W_D1); 
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * D1_coupling(i, j);
        };
        
        // Options or evaluating only one piece
        static const int kSwaveOnly = 1;
        static const int kDwaveOnly = 2;

        // -----------------------------------------------------------------------
        private:

        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            double gs = (_option == kDwaveOnly) ? 0. : H1_S;
            double gd = (_option == kSwaveOnly) ? 0. : H1_D;

            double swave = gs*sqrt(_mpc*_mpc + _mc2)            * delta(i,j);
            double dwave = gd*_mpc*_mpc*(3.*phat_1(i)*phat_1(j) - delta(i,j));

            return sqrt(M_D1*M_DSTAR)*(swave + dwave);
        };

        inline void recalculate()
        {
            // Update arguments with floating Y and Z meson masses for the triangle
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _AD = _T.eval();

            // Gather all couplings
            _AD *= _Y->coupling()/sqrt(2.) * sqrt(M_D1*M_D*_Y->mass());  
            // Skip D-wave D1 -> D* pi which we factor out
            _AD *= _Zc->coupling()         * sqrt(M_D*M_DSTAR*_Zc->mass()); 
            _AD *= _Zc->propagator(_sab);
            _AD *= _Zc->coupling()         * sqrt(M_D*M_DSTAR*_Zc->mass()); 
        };

        // Energy dependent D wave strength
        complex _AD;  

        // Scalar triangle function
        hadMolee::triangle _T;

        // This channel recieves contribution from the Z(3900) & Y (molecular) couplings
        molecule _Zc, _Y;
    };

    // ---------------------------------------------------------------------------
    // Also include the D1 tree transition 
    // No free parameters since this assume to be known background

    class tree : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        tree(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_tree")
        : amplitude_base(key, xkinem, V, 0, "DsDpi_dwave", id),
          _D1(breit_wigner::kNonrelativistic, M_D1, W_D1),
          _Y(get_molecular_component(_V))
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * D1_coupling(i, j);
        };

        // Options or evaluating only one piece
        static const int kDefault   = 0;
        static const int kSwaveOnly = 1;
        static const int kDwaveOnly = 2;

        // -----------------------------------------------------------------------
        private:

        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            double gs = (_option == kDwaveOnly) ? 0. : H1_S;
            double gd = (_option == kSwaveOnly) ? 0. : H1_D;

            double swave = gs*sqrt(M_PION*M_PION + _mpc*_mpc)   * delta(i,j);
            double dwave = gd*_mpc*_mpc*(3.*phat_1(i)*phat_1(j) - delta(i,j));

            return sqrt(M_D1*M_DSTAR)*(swave + dwave);
        };

        inline void recalculate()
        {            
            _AD  = _Y->coupling() / sqrt(2.) * sqrt(_V->mass()*M_D1*M_D);  
            _AD *= _D1.eval(_sac); // D1 propagator
            // Skip D1 -> D* pi which we factor out and put above
        };

        // Couplings
        complex _AD;  // Energy dependent D wave strength

        // In addition we have the tree level transition
        breit_wigner _D1;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
    };
};

#endif