// Specific implementations of amplitudes relevant for the D* D pi final state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DSDPI_HPP
#define DSDPI_HPP

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
    // Contact-like interaction in S-wave
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients
    class swave : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        swave(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "DsDpi_swave")
        : amplitude_base(key, xkinem, Y, 2, "DsDpi_swave", id),
          _Zc(make_molecule<DsD_molecule>())
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * delta(i,j);
        };

        // -----------------------------------------------------------------------
        private:

        // Parameterize with two real polynomial coefficients
        // Values from Leon's note
        double _a = 2.501, _b = -15.12;

        // S-wave coupling strength
        complex _AS; 

        // This channel recieves contribution from the Z(3900)
        molecule _Zc;

        // S-wave amplitude is easy, because its only the Z propagator that needs to be calcualted
        inline void recalculate()
        {
            // Multiply by the helicity frame energy to make sure amplitude respects Goldstone theorem.
            double omega_pi = sqrt(M_PION*M_PION + _mpc*_mpc);

            // sab is assumed to be DsD channel
            _AS  = _a*(_b + _sab)*_Zc->propagator(_sab)*omega_pi;
        };
    };

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

        // D1 -> D* pi coupling
        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            return sqrt(M_D1*M_DSTAR)*_mpc*_mpc*(H1_D*(3.*phat_1(i)*phat_1(j) - delta(i,j)) + H1_S*delta(i,j));
        };
        
        inline void recalculate()
        {
            // Y mass and couplings
            double y   = _Y->molecular_coupling();
            double M_Y = _V->pole_mass();

            // Update arguments with floating Y and Z meson masses for the triangle
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _AD = _T.eval();

            // This gets multiplied by the propagator of the Z prop and triangle function
            double  z   = _Zc->molecular_coupling();
            double  M_Z = M_ZC3900;
            complex G_Z = _Zc->propagator(_sab);

            // Couplings at the vertices of the triangle
            _AD *= y/sqrt(2.) * sqrt(M_Y*M_D1*M_D);  
            // Skip D-wave D1 -> D* pi which we factor out
            _AD *= z          * sqrt(M_D*M_DSTAR*M_Z); 

            // Z decay vertex
            _AD *= z * G_Z    * sqrt(M_D*M_DSTAR*M_Z);
        };

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
            return sqrt(M_D1*M_DSTAR)*_mpc*_mpc*(H1_D*(3.*phat_1(i)*phat_1(j) - delta(i,j)) + H1_S*delta(i,j));
        };
        
        inline void recalculate()
        {            
            // Y mass and couplings
            double y   = _Y->molecular_coupling();
            double M_Y = _V->pole_mass();

            // D1 propagator
            complex G_D1 = _D1.eval(sqrt(_sac));

            // Only two verices
            _AD  = G_D1; 
            _AD *= y/sqrt(2.) * sqrt(M_Y*M_D1*M_D);  
            // Skip D-wave D1 -> D* pi which we factor out
        };

        // Couplings
        complex _AD;  // Energy dependent D wave strength

        // In addition we have the tree level transition
        breit_wigner _D1;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
    };

    // ---------------------------------------------------------------------------
    // Contact-like interaction involving the Psi(4160)
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients
    class psi_contact : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        psi_contact(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "DsDpi_psi")
        : amplitude_base(key, xkinem, Y, 0, "DsDpi_psi", id),
          _Zc(make_molecule<DsD_molecule>())
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * delta(i,j);
        };

        // -----------------------------------------------------------------------
        private:

        double _mod_psi  = 1.9;

        // S-wave coupling strength
        complex _AS; 

        // This channel recieves contribution from the Z(3900)
        molecule _Zc;

        // S-wave amplitude is easy, because its only the Z propagator that needs to be calcualted
        inline void recalculate()
        {
            // Multiply by the helicity frame energy to make sure amplitude respects Goldstone theorem.
            double omega_pi = sqrt(M_PION*M_PION + _mpc*_mpc);

            // Add the psi(4160) part by removing the Ylineshape
            _AS = _mod_psi * _Zc->propagator(_sab)*omega_pi;
        };
    };
};

#endif