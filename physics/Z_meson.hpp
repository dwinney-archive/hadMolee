// Model for the hadronic molecule content of the D*D meson relevant for the Zc(3900)
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef Z_MESON_HPP
#define Z_MESON_HPP

#include "molecule.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // D* D molecule relevant for the Z meson
    // this class is only molecular because its axial-vector and doesnt couple to the photon

    class DsD_molecule : public molecular
    {
        // -----------------------------------------------------------------------
        public:

        DsD_molecule()
        : molecular(M_DSTAR, M_D)
        {
            // Given a pole mass calculate the bare mass
            _pole_mass           = 3.87;
            _z                   = 2.871;
            _molecular_coupling  = _z;

            _sigma_pole          = _z*_z*self_energy(_pole_mass*_pole_mass);
            _bare_mass           = _pole_mass - real(_sigma_pole);

            // Calcualte the non-moleculat component to the width from the total
            _total_width         = 100E-3;
            _nonmol_width        = _total_width - imag(_sigma_pole);
        };

        // The propagator gains contributions from the self-energy
        inline complex propagator(double s)
        {
            double E = sqrt(s);
            complex D = E - _bare_mass - _z*_z*self_energy(s) + I*_nonmol_width/2;
            
            return I / D;
        };  

        private:
        
        // Rename molecular coupling to z
        double _z;
    };
};

#endif