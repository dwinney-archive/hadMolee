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
            // Given a pole mass calculate the bare mass
            _pole_mass           = 4.23;
            _y                   = 42.8E-3;
            _molecular_coupling  = _y;

            _sigma_pole          = _y*_y*self_energy(_pole_mass*_pole_mass);
            _bare_mass           = _pole_mass + real(_sigma_pole);

            // Calcualte the non-moleculat component to the width from the total
            _total_width         = 46E-3;
            _nonmol_width        = _total_width + 2*imag(_sigma_pole);
        };

        // The propagator gains contributions from the self-energy
        complex propagator(double s)
        {
            double  E = sqrt(s);
            complex D = E - _pole_mass - _y*_y*self_energy(s) + I*_nonmol_width/2.;
            return I / (2*D);
        };

        // Since Y-meson is also charmonium-like it requires a photon coupling
        complex photon_coupling()
        {
            return I * E * _pole_mass*_pole_mass / _fY;
        };

        private:

        // Rename the molecular coupling to just y
        double _y;

        // Y decay constant
        double _fY = 2.194;
    };
};

#endif