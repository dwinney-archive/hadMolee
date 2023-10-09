// Conventional charmonium lineshape model 
//
// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef CHARMONIUM_HPP
#define CHARMONIUM_HPP

#include "lineshape.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "constants.hpp"

namespace hadMolee
{
    class charmonium : public charmoniumlike
    {
        public:

        charmonium(std::vector<double> pars, std::string id = "psi(4160)")
        : charmoniumlike(3, id)
        {
            _mass  = pars[0];
            _width = pars[1]; 
            _fpsi  = decay_constant(pars[2]);
        };

        inline complex propagator(double s)
        {
            return 1./(s - _mass*_mass + I*_mass*_width);
        };

        inline complex photon_coupling()
        {
            return E * _mass*_mass / _fpsi;
        };

        inline double pole_mass(){ return _mass; };

        private:

        inline double decay_constant(double BR)
        {
            double Gamma_ee = _width * BR;
            return sqrt( 4.*PI*ALPHA* ALPHA * _mass / (3. * Gamma_ee));
        };

        double _fpsi  = 0., _mass = 0., _width = 0.;
    };
};

#endif