// Model for the hadronic molecule content of the D*D meson relevant for the Zc(3900)
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef Y_MESON_HPP
#define Y_MESON_HPP

#include "molecule.hpp"
#include "lineshape.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Implementation of D1 D molecule for the Y(4260)
    // We inherit from both charmoniumlike and molecular to define a lineshape

    class Y_meson : public charmonium, public molecular
    {
        // -----------------------------------------------------------------------
        public:

        Y_meson()
        : charmonium(4), molecular(M_D1, M_D)
        {};

        // Setting parameters involves both allocating lineshape and molecule pars
        inline void set_parameters(std::vector<double> pars)
        {
            charmonium::check_size(pars);
            charmonium::_mass    = pars[0];
            // Dont need charmonium::_width 
            charmonium::_fV      = pars[3];
            
            molecular::_mass     = pars[0];
            molecular::_nm_width = pars[1];
            molecular::_coupling = pars[2];
            
            // When updating mass also update self-energy subtraction
            _reSigmaPole = std::real(self_energy(pars[0]*pars[0]));
        };

        // Override charmonium::propagator to always use the molecular one
        inline complex propagator(double s){ return molecular::propagator(s); };
    };
};

#endif