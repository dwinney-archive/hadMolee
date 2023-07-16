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
#include "Y_meson.hpp"
#include "Z_meson.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Contact-like interaction in S-wave
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients
 class DsDpi_swave : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        DsDpi_swave(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "DsDpi_swave")
        : amplitude_base(key, xkinem, Y, 2, "DsDpi_swave", id)
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * delta(i,j);
        };

        // Set the two parameters
        // These are the coefficients of the linear polynomial multiplying the s-wave strength
        inline void set_parameters(std::vector<double> par)
        {
            check_nParams(par);

            _a = par[0];
            _b = par[1];
        };

        // -----------------------------------------------------------------------
        private:

        // Parameterize with two real polynomial coefficients
        // Values from Leon's note
        double _a = 2.501, _b = -15.12;

        // S-wave coupling strength
        complex _AS; 

        // This channel recieves contribution from the Z(3900)
        molecule _Zc = make_molecule<DsD_molecule>();

        // S-wave amplitude is easy, because its only the Z propagator that needs to be calcualted
        inline void recalculate()
        {
            // Multiply by the helicity frame energy to make sure amplitude respects Goldstone theorem.
            double p_pi     = _kinematics->decay_momentum_c(_s, _sab);
            double omega_pi = sqrt(M_PION*M_PION + p_pi*p_pi);

            // sab is assumed to be DsD channel
            _AS = _a * (_b + _sab) * _Zc->propagator(_sab) * omega_pi;
        };
    };

    // ---------------------------------------------------------------------------
    // Tree level transition through intermediate D1
    // No free parameters since this assume to be known background

    class DsDpi_tree : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        DsDpi_tree(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_tree")
        : amplitude_base(key, xkinem, V, 0, "DsDpi_tree", id)
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return 1E3 * _AD * (3.*phat(i)*phat(j) - delta(i,j));
        };

        // -----------------------------------------------------------------------
        private:

        inline void recalculate()
        {
            double p_pion  = _kinematics->decay_momentum_c(_s, _sab);
            double p2_pion = p_pion*p_pion;

            // sab is assumed to be DsD channel
            _AD = - sqrt(2./3.)*_hD/_fpi * p2_pion;

            // Multiply by the propagator of the D1 and y coupling
            _AD *= _Y->molecular_coupling() / sqrt(2) /* * _D1.eval(_sac) */;
        };
        
        // Energy dependent D wave strength
        complex _AD;                         

        // HQSS couplings
        double  _fpi = sqrt(2) * 92E-3;    // Pion decay constant in GeV
        double  _hD  = 0.89;               // HQSS constant in GeV-1

        // This channel has a D1 resonance in the Ds pi channel
        breit_wigner _D1 = breit_wigner(breit_wigner::kRelativistic, M_D1, W_D1);

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y = get_molecular_component(_V);
    };

    // // ---------------------------------------------------------------------------
    // // Transition involving D1DsD triangle and Z meson
    // // No free parameters since this assume to be known background

    // class DsDpi_triangle : public amplitude_base
    // {
    //     // -----------------------------------------------------------------------
    //     public:
        
    //     // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
    //     // We default to the nonrel version
    //     DsDpi_triangle(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_triangle")
    //     : amplitude_base(key, xkinem, V, 0, "DsDpi_triangle", id)
    //     {};

    //     // The reduced amplitude corresponds to the S-wave contact-like interaction
    //     // however it receives contributions from the propagator of the Z meson
    //     inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
    //     {
    //         // Being D-wave we get the appropriate projector
    //         // Assume pion (particle c) defined the +z direction
    //         return _AS * s_wave(i,j) + _AD * d_wave(i,j);
    //     };

    //     // -----------------------------------------------------------------------
    //     private:
        
    //     inline void recalculate()
    //     {
    //         double p_pion = _kinematics->decay_momentum_c(_s, _sab);

    //         // sab is assumed to be DsD channel
    //         _AS = I * _hS * sqrt(p_pion*p_pion + _mc2);
    //         _AD = I * _hD * p_pion*p_pion;

    //         // Update arguments with floating Y and Z meson masses
    //         _T.set_internal_masses(_internal); 
    //         _T.set_external_masses({_W, sqrt(_sab), M_PION});
    //         _T.add_width(b, W_D1); 
            
    //         // Multiply by the propagator of the Z and triangle function
    //         double z = _Zc->molecular_coupling();
            
    //         // The normalization is matched to Qiang's paper
    //         complex T = _internal[0] * _internal[1] * _internal[2] * _T.eval();

    //         // Common pieces for both S and D wave
    //         _AS *= - (_Y->molecular_coupling() / sqrt(2.)) * z*z * T * _Zc->propagator(_sab);
    //         _AD *= - (_Y->molecular_coupling() / sqrt(2.)) * z*z * T * _Zc->propagator(_sab);
    //     };

    //     // Couplings
    //     double _hS = 0., _hD = HPRIME_UPPER / F_PION;  // D1 -> D*pi coupling for the S-wave and the D-wave
    //     complex _AS, _AD;                              // Energy dependent S and D wave strengths

    //     // On-shell masses involved in the triangle
    //     std::array<double,3> _internal = {M_DSTAR, M_D1, M_D};

    //     // Default to using nonrelativistic version of the triangle
    //     triangle _T = triangle(nonrelativistic);

    //     // This channel recieves contribution from the Z(3900)
    //     molecule _Zc = make_molecule<DsD_molecule>();

    //     // It also explciity depends on D1D molecular nature of the Y state
    //     molecule _Y = get_molecular_component(_V);
    // };
};

#endif