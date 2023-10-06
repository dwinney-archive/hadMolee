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

        complex propagator(double s)
        {
            complex D = sqrt(s) - _mass + I*_width/2.;
            return 1. / (2. * _mass * D);
        };

        complex photon_coupling()
        {
            return E * _mass*_mass / _fpsi;
        };

        double pole_mass(){ return _mass; };

        private:

        virtual double decay_constant(double BR)
        {
            return 1.;
        };

        double _fpsi  = 0., _mass = 0., _width = 0.;
    };
};

#endif