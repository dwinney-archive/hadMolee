// Conventional charmonium lineshape model for the psi(4160) 
//
// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef PSI4160_HPP
#define PSI4160_HPP

#include "lineshape.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "Z(3900).hpp"
#include "constants.hpp"

namespace hadMolee
{
    class psi4160 : public charmoniumlike
    {
        public:

        psi4160(std::string id = "psi(4160)")
        : charmoniumlike(0, id)
        {};

        virtual complex propagator(double s)
        {
            complex D = sqrt(s) - _mass + I*_width/2.;
            return 1. / (2. * _mass * D);
        };

        virtual complex photon_coupling()
        {
            return E * _mass*_mass / _fpsi;
        };

        private:

        double _fpsi  = 43.8427;
        double _mass  = 4.191;
        double _width = 70.E-3;
    };

    // ---------------------------------------------------------------------------
    // Contact-like interaction in S-wave
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients
    class DsDpi_psi : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        DsDpi_psi(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "DsDpi_psi")
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
            double p_pi     = _kinematics->decay_momentum_c(_s, _sab);
            double omega_pi = sqrt(M_PION*M_PION + p_pi*p_pi);

            // sab is assumed to be DsD channel

            // Add the psi(4160) part by removing the Ylineshape
            _AS = _mod_psi * _Zc->propagator(_sab)*omega_pi;
        };
    };
};

#endif