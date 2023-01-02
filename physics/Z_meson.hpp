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

    class DsD_molecule : public molecular
    {
        // -----------------------------------------------------------------------
        public:

        DsD_molecule()
        : molecular(M_DSTAR, M_D)
        {
            // Mass and Width from PDG
            _pole_mass          = M_ZC3900;
            _total_width        = W_ZC3900;
            
            // Coupling taken from [1]
            _molecular_coupling      = ZBARE_QQ2016;  
            
            // residual width taken to recover the full PDG width at the pole
            _nonmol_width = _total_width - 2.* imag(self_energy(_pole_mass));
        };

        // The propagator gains contributions from the self-energy
        inline complex propagator(double s)
        {
            double E = sqrt(s);
            double z = _molecular_coupling;

            complex D = E - _pole_mass + I * (z*z*self_energy(E) + _nonmol_width/2.);
            
            return I / (2.*D);
        };  

        // Self-energy from bubble diagram of D* D scattering 
        inline complex self_energy(double E)
        {
            double eps = mass_difference(E);
            double mu  = reduced_mass();
            
            return (1. / (8.*PI)) * csqrt(2.*mu*mu*mu*std::abs(eps)) * ( (eps>=0) + I*(eps<0) );
        };

        // -----------------------------------------------------------------------
        private:

        // Total width of the Z from the PDG
        double _total_width;
    };
};

#endif