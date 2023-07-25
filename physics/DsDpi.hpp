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
        : amplitude_base(key, xkinem, Y, 2, "DsDpi_swave", id),
          _Zc(make_molecule<DsD_molecule>())
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * sqrt(sqrt(_s)*M_DSTAR*M_D) * delta(i,j);
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
            double p_pi     = _kinematics->decay_momentum_c(_s, _sab);
            double omega_pi = sqrt(M_PION*M_PION + p_pi*p_pi);

            // sab is assumed to be DsD channel
            _AS = _a*(_b + _sab)*_Zc->propagator(_sab)*omega_pi;
        };
    };

    // ---------------------------------------------------------------------------
    // Transition involving D1DsD triangle with the Z as well as 
    // No free parameters since this assume to be known background

    class DsDpi_triangle : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        DsDpi_triangle(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_triangle")
        : amplitude_base(key, xkinem, V, 0, "DsDpi_dwave", id),
          _T(triangle::kNonrelativistic), 
          _Zc(make_molecule<DsD_molecule>()),
          _Y(get_molecular_component(_V))
        {
            _T.set_internal_masses(_internal); 
            _T.add_width(2, W_D1); 
        };

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * sqrt(sqrt(_s)*M_DSTAR*M_D) *_ppi*_ppi*(3.*phat(i)*phat(j) - delta(i,j));
        };

        // -----------------------------------------------------------------------
        private:
        
        inline void recalculate()
        {
            // Update arguments with floating Y and Z meson masses for the triangle
            _T.set_external_masses({_W, sqrt(_sab), M_PION});

            // Update the pion momentum 
            _ppi = _kinematics->decay_momentum_c(_s, _sab);
            
            // This gets multiplied by the propagator of the Z prop and triangle function
            double z = _Zc->molecular_coupling();

            _AD  = I * sqrt(2./3.)*_hp/_fpi;
            _AD *= (_Y->molecular_coupling()/sqrt(2.));
            _AD *= z*z*M_D1*M_D*M_DSTAR*(2*_T.eval())*_Zc->propagator(_sab);
            // Factor of 2 comes from sum of charge conjugated diagrams
        };

        // Couplings
        complex _AD;  // Energy dependent D wave strength

        // Couplings related to pion
        double  _ppi;                            // Pion 3-momentum
        double  _fpi = /* sqrt(2)* */92.8E-3;    // Pion decay constant in GeV
        double  _hp  =  /* 0.89 */ 0.62;         // HQSS constant in GeV-1   

        // On-shell masses involved in the triangle
        std::array<double,3> _internal = {M_DSTAR, M_D1, M_D};

        // Default to using nonrelativistic version of the triangle
        triangle _T;

        // This channel recieves contribution from the Z(3900)
        molecule _Zc;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
    };

    // ---------------------------------------------------------------------------
    // Also include the D1 tree transition 
    // No free parameters since this assume to be known background

    class DsDpi_tree : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        DsDpi_tree(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "DsDpi_tree")
        : amplitude_base(key, xkinem, V, 0, "DsDpi_dwave", id),
          _D1(breit_wigner::kNonrelativistic,  M_D1, W_D1),
          _Y(get_molecular_component(_V))
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _AD * sqrt(sqrt(_s)*M_DSTAR*M_D) *_ppi*_ppi*(3.*phat(i)*phat(j) - delta(i,j));
        };

        // -----------------------------------------------------------------------
        private:
        
        inline void recalculate()
        {
            // Update the pion momentum 
            _ppi = _kinematics->decay_momentum_c(_s, _sab);
            
            _AD  = sqrt(2./3.)*_hp/_fpi;
            _AD *= (_Y->molecular_coupling()/sqrt(2.));
            _AD *= M_D1*_D1.eval(sqrt(_sac));
        };

        // Couplings
        complex _AD;  // Energy dependent D wave strength

        // Couplings related to pion
        double  _ppi;                            // Pion 3-momentum
        double  _fpi = /* sqrt(2)* */92.8E-3;    // Pion decay constant in GeV
        double  _hp  =  /* 0.89 */ 0.62;         // HQSS constant in GeV-1   

        // In addition we have the tree level transition
        breit_wigner _D1;

        // It also explciity depends on D1D molecular nature of the Y state
        molecule _Y;
    };
};

#endif