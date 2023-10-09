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

    class D1D_molecule : public charmoniumlike, public molecular
    {
        // -----------------------------------------------------------------------
        public:

        D1D_molecule(std::string id = "Y(4260)")
        : charmoniumlike(4, id), molecular(M_D1, M_D)
        {
            _pole_mass           = 4.2281;
            _y                   = 0.1145 / sqrt(_pole_mass * M_D1 * M_DSTAR);
            _total_width         = 42.6E-3;
            _fY                  = 0.672;
            _molecular_coupling  = _y;

            _bare_mass           = _pole_mass;
            _nonmol_width        = _total_width;

            _reSigmaPole = std::real(self_energy(_bare_mass*_bare_mass));
        };

        inline double pole_mass(){ return _pole_mass; };

        // The propagator gains contributions from the self-energy
        inline complex propagator(double s)
        {
            complex Sigma = self_energy(s) - _reSigmaPole;
            complex D = sqrt(s) - _pole_mass - _y*_y*Sigma + I*_nonmol_width/2.;
            
            return 1. / (2. * _pole_mass * D);
        };  

        // Since Y-meson is also charmonium-like it requires a photon coupling
        complex photon_coupling()
        {
            return E * _pole_mass*_pole_mass / _fY;
        };

        private:

        double _y, _fY, _reSigmaPole;
    };
};

#endif