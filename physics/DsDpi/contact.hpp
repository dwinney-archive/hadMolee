// Specific implementations of amplitudes relevant for the D* D pi final state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DSDPI_CONTACT_HPP
#define DSDPI_CONTACT_HPP

#include <memory>

#include "kinematics.hpp"
#include "breit_wigner.hpp"
#include "amplitude.hpp"
#include "molecule.hpp"
#include "constants.hpp"
#include "triangle.hpp"

namespace hadMolee::DsDpi
{
      // ---------------------------------------------------------------------------
    // Contact-like interaction in S-wave
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients
    class contact : public amplitude_base
    {
        // -----------------------------------------------------------------------
        public:
        
        contact(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "DsDpi_psi")
        : amplitude_base(key, xkinem, Y, 2, "DsDpi_psi", id),
          _Zc(make_molecule(M_D, M_DSTAR))
        {
            _Zc->set_parameters({3.9, 50.E-3, 4.66/sqrt(M_DSTAR*M_D*3.9)});
        };

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            complex decay_amplitude =  _AS * sqrt(M_PION*M_PION + _mpc*_mpc);
            
            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * sqrt(M_PION*M_PION + _mpc*_mpc) * delta(i,j);
        };

        inline void set_parameters( std::vector<double> pars )
        {
            check_nParams(pars);
            _a = pars[0];
            _b = pars[1];
        };

        // -----------------------------------------------------------------------
        private:

        // Amplitude constants
        double _a  = 0., _b = 0.;

        // S-wave coupling strength
        complex _AS; 

        // This channel recieves contribution from the Z(3900)
        molecule _Zc;

        // S-wave amplitude is easy, because its only the Z propagator that needs to be calcualted
        inline void recalculate()
        {
            _AS  = _a*(_b + _sab) * _Zc->propagator(_sab);
        };
    };
};

#endif