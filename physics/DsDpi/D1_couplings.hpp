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
#include "lineshapes/Y(4260).hpp"
#include "lineshapes/Z(3900).hpp"

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
          _Zc(make_molecule<DsD_molecule>()),
          _Y(get_molecular_component(_V))
        {
            _T.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _T.add_width(2, W_D1); 
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * D1_coupling(i, j);
        };

        // -----------------------------------------------------------------------
        private:

        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            double gs = (_debug == 1) ? 0. : H1_S;
            double gd = (_debug == 2) ? 0. : H1_D;

            double swave = gs*_omega_pi*delta(i,j);
            double dwave = gd*_mpc*_mpc*(3.*phat_1(i)*phat_1(j) - delta(i,j));

            return sqrt(M_D1*M_DSTAR)*(swave + dwave);
        };

        inline void recalculate()
        {
            // Y mass and couplings
            double gy   = _Y->molecular_coupling();
            double M_Y  = _V->pole_mass();

            // Update arguments with floating Y and Z meson masses for the triangle
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _AD = _T.eval();

            // Pion energy for S-wave
            _omega_pi = sqrt(_mpc*_mpc + _mc2);

            // This gets multiplied by the propagator of the Z prop and triangle function
            double  gz   = _Zc->molecular_coupling();
            double  M_Z = M_ZC3900;
            complex G_Z = _Zc->propagator(_sab);

            // Couplings at the vertices of the triangle
            _AD *= gy/sqrt(2.) * sqrt(M_Y*M_D1*M_D);  
            // Skip D-wave D1 -> D* pi which we factor out
            _AD *= gz          * sqrt(M_D*M_DSTAR*M_Z); 

            // Z decay vertex
            _AD *= gz * G_Z    * sqrt(M_D*M_DSTAR*M_Z);
        };

        // Pion energy
        double _omega_pi;

        // Energy dependent D wave strength
        complex _AD;  

        // Scalar triangle function
        hadMolee::triangle _T;

        // This channel recieves contribution from the Z(3900)
        molecule _Zc;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
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
          _D1(breit_wigner::kNonrelativistic,  M_D1, W_D1),
          _Y(get_molecular_component(_V))
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * D1_coupling(i, j);
        };

        // -----------------------------------------------------------------------
        private:

        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            double gs = (_debug == 1) ? 0. : H1_S;
            double gd = (_debug == 2) ? 0. : H1_D;

            double swave = gs*_omega_pi*delta(i,j);
            double dwave = gd*_mpc*_mpc*(3.*phat_1(i)*phat_1(j) - delta(i,j));

            return sqrt(M_D1*M_DSTAR)*(swave + dwave);
        };

        inline void recalculate()
        {            
            // Y mass and couplings
            double y   = _Y->molecular_coupling();
            double M_Y = _V->pole_mass();

            // D1 propagator
            complex G_D1 = _D1.eval(sqrt(_sac));
            
            // Pion energy for S-wave
            _omega_pi = sqrt(_mpc*_mpc + _mc2);

            // Only two verices
            _AD  = G_D1; 
            _AD *= y/sqrt(2.) * sqrt(M_Y*M_D1*M_D);  
            // Skip D-wave D1 -> D* pi which we factor out
        };

        // Couplings
        complex _AD;  // Energy dependent D wave strength

        // Pion energy
        double _omega_pi;

        // In addition we have the tree level transition
        breit_wigner _D1;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
    };
};

#endif