// Specific implementations of amplitudes relevant for the D* D pi final state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DSDPI_SWAVE_HPP
#define DSDPI_SWAVE_HPP

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
};

#endif