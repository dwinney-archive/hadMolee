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
            _pole_mass           = 3.9;
            _z                   = 4.66 / sqrt(M_DSTAR * M_D * _pole_mass);
            _total_width         = 50.E-3;
            _molecular_coupling  = _z;

            // IF WE TAKE THE INPUT MASS AS THE BARE MASS
            _bare_mass    = _pole_mass;
            _nonmol_width = _total_width;

            _reSigmaPole = std::real(self_energy(_bare_mass*_bare_mass));
        };

        // The propagator gains contributions from the self-energy
        inline complex propagator(double s)
        {
            double  E     = sqrt(s);
            complex Sigma = self_energy(s) - _reSigmaPole;
            complex D = E - _pole_mass - _z*_z*Sigma + I*_nonmol_width/2.;
            
            return 1. / (2. * _pole_mass * D);
        };  

        private:
        
        // Rename molecular coupling to z
        double _z;
        double _reSigmaPole;
    };
};

#endif