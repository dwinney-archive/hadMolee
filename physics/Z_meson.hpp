// Model for the hadronic molecule content of the D*D meson relevant for the Zc(3900)
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef Z_MESON_HPP
#define Z_MESON_HPP

#include "hadronic_molecule.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // D* D molecule relevant for the Z meson

    class DsD_molecule : public hadronic_molecule
    {
        // -----------------------------------------------------------------------
        public:

        DsD_molecule()
        : hadronic_molecule(M_DSTAR, M_D)
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
        inline std::complex<double> propagator(double s)
        {
            double E = sqrt(s);
            double z = _molecular_coupling;

            std::complex<double> D = E - _pole_mass + XI * (z*z*self_energy(E) + _nonmol_width/2.);
            
            return XI / (2.*D);
        };  

        // Self-energy from bubble diagram of D* D scattering 
        inline std::complex<double> self_energy(double E)
        {
            double eps = mass_difference(E);
            double mu  = reduced_mass();
            
            return (1. / (8.*PI)) * sqrt(2.*mu*mu*mu*std::abs(eps)) * ( XR*(eps>=0) + XI*(eps<0) );
        };

        // -----------------------------------------------------------------------
        private:

        // Total width of the Z from the PDG
        double _total_width;
    };
};

#endif