// Specific implementations of amplitudes relevant for the Jpsi pi pi final state

// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef JPSIPIPI_COUNTERTERMS_HPP
#define JPSIPIPI_COUNTERTERMS_HPP

#include <memory>

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "triangle.hpp"
#include "Y(4260).hpp"

namespace hadMolee::Jpsipipi
{
    class triangle_CT : public amplitude_base
    {
        public:
        
        triangle_CT(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "triangle_CT")
        : amplitude_base(key, xkinem, Y, 0, "triangle_CT", id), 
        _Y(get_molecular_component(Y)), _Zc(make_molecule(M_DSTAR, M_D)),
        _T(triangle::kLoopTools)
        {
            // Set up Zc propagator
            _Zc->set_parameters({3.9, 50.E-3, 4.66/sqrt(M_DSTAR*M_D*3.9)});
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index k)
        {
            complex result = 0.;
            for (auto j : C_INDICES)
            {
                // particle c couples to the D1 vertex
                _pi1 = particle::c;
                result += D1_coupling(i, j) * T1() * delta(j, k);

                // particle b couples to the D1 vertex
                _pi1 = particle::b;
                result -= D1_coupling(i, j) * T1() * delta(j, k);
            };
            return  _C * result / sqrt(2.);
        };

        protected:

        // Save which particle is considered pi1, by default this is particle c
        particle _pi1 = particle::c;

        // D1 -> D* pi coupling
        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            return sqrt(M_D1*M_DSTAR)*(H1_D*modp()*modp()*(3.*phat_1(i)*phat_1(j) - delta(i,j)) + H1_S*delta(i,j));
        };

        inline void recalculate()
        {
            // Y mass and couplings
            double gy  = _Y->coupling();
            double M_Y = sqrt(_s);
            double gz  = _Zc->coupling();
            double M_Z = M_ZC3900;

            // Update floating masses in triangle 1 and evaluate
            // Also include the Zc propagator and omega_pi here
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _T.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _T.add_width(2, W_D1); 
            _T1[0] = _T.eval() * _Zc->propagator(_sab) * sqrt(_mpb*_mpb + _mb2);
            // Swap pions 
            _T.set_external_masses({_W, sqrt(_sac), M_PION});
            _T1[1] = _T.eval() * _Zc->propagator(_sac) * sqrt(_mpc*_mpc + _mc2);

            // Couplings at the vertices of the first triangle
            _C  = gy/sqrt(2.) * sqrt(M_Y*M_D1*M_D);    // Y -> D1 D
            // Skip the D1 -> D* pi coupling which we include above
            _C *= gz          * sqrt(M_D*M_DSTAR*M_Z); // D* D -> Z 
            _C *= CT_TRI; 
        };

        // Aliases for the pion momenta specifying which is pi1 and pi2
        inline double modp(){ return (_pi1 == particle::c) ? _mpc : _mpb; };

        // Returns T1 for current permutation of pions
        inline complex T1(){  return (_pi1 == particle::c) ? _T1[0] : _T1[1]; };

        // Scalar prefactors 
        double _C     = 0.;

        // Triangle functions, this is used to evaluate both triangles
        triangle _T;

        // Cached quantities related to the triangle
        // index 0 couples particle c to the D-wave vertex while 1 is particle b
        std::array<complex,2> _T1;  // Scalar Y -> Zc pi triangle

        // We need access the molecular nature of the Y and Zc state
        molecule _Y, _Zc;
    };
};

#endif