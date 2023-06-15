// Model for Y meson, this includes both the molecular content as a D1D state
// as well as its coupling to the photon being a 1-- state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef Y_MESON_HPP
#define Y_MESON_HPP

#include "lineshape.hpp"
#include "molecule.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Implementation of D1 D molecule for the Y(4260)

    class D1D_molecule : public charmoniumlike, public molecular
    {
        // -----------------------------------------------------------------------
        public:

        D1D_molecule(std::string id = "Y(4260)")
        : charmoniumlike(4, id), molecular(M_D1, M_D)
        {
            // Set up the derivator 
            wsigma = ROOT::Math::Functor1D(this, &D1D_molecule::resigma);
            dsigma.SetFunction(wsigma);

            // Default parameters use the QQ2016 values
            set_parameters({MY_QQ2016, YBARE_QQ2016, YNM_WIDTH_QQ2016, F_Y_QQ2016});
        };

        // The propagator gains contributions from the self-energy
        complex propagator(double s)
        {
            double E = sqrt(s);
            complex D = E - _pole_mass - _Z*self_energy(E) + I*_nonmol_width/2.;
            return I * _Z / (2. * D);
        };

        // Since Y-meson is also charmonium-like it requires a photon coupling
        complex photon_coupling()
        {
            return I * E * _pole_mass*_pole_mass / _fY;
        };

        // When we change the pole mass, we must recalculate renormalization quantities
        // Set pole mass ,coupling, and non-mol width in a single call
        inline void set_parameters(std::vector<double> pars)
        {
            check_size(pars);

            _pole_mass              = pars[0];
            _molecular_coupling     = pars[1];
            _nonmol_width           = pars[2];
            _fY                     = pars[3];

            // Precalculate relevant quantities
            _reS  = resigma(_pole_mass);
            _redS = dsigma.Eval(_pole_mass);
            _Z    = 1. / (1. - _redS);
        };

        private:

        // // Self-energy from bubble of D1 D scattering and dressed with elastic scattering
        // // renomalized
        // complex self_energy(double E)
        // {
        //     return sigma(E) - _reS - (E - _pole_mass) * _redS;
        // };

        // Bare self-energy just from the bubble of D1 D scattering
        inline complex sigma(double E)
        {
            double eps = mass_difference(E);
            double mu  = reduced_mass();
            double y   = _molecular_coupling;

            return (-I/ (8.*PI)) * csqrt( 2.*mu*mu*mu* (eps + I*W_D1/2.) ) * (y*y);
        };
        inline double resigma(double E){ return real(sigma(E)); };

        // Need to be able to calculate the derivative of the above self-energy
        ROOT::Math::Functor1D wsigma;
        ROOT::Math::Derivator dsigma;

        double _Z, _reS, _redS, _fY;
    };
};

#endif