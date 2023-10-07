// Specific implementations of amplitudes relevant for the Jpsi pi pi final state

// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef JPSIPIPI_BOXES_HPP
#define JPSIPIPI_BOXES_HPP

#include <memory>

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "box.hpp"
#include "lineshapes/Y(4260).hpp"

namespace hadMolee::Jpsipipi
{
    class Box_I : public amplitude_base
    {
        public:

        Box_I(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box I")
        : amplitude_base(key, xkinem, Y, 0, "Jpsipipi_BoxI", id), 
          _B(box::kLoopTools), _Y(get_molecular_component(Y))
        {
            _B.set_internal_masses(_internal);
            _B.set_external_masses(_external);
            _B.add_width(0, W_D1);
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index k)
        {
            complex result = 0.;
            for (auto j : C_INDICES)
            {
                result += (H1_D*_mpc*_mpc*(3.*p_c(i)*p_c(j) - delta(i,j)) + H1_S*delta(i,j)) * M(j, k);
            };
            return  result;
        };

        private:

        inline void recalculate()
        {
            // Y mass and couplings
            double gy   = _Y->molecular_coupling();
            double M_Y = sqrt(_s);

            // Modulus of pion 3-momenta 
            _mpc = _kinematics->decay_momentum_c(_s, _sab);
            _mpb = _kinematics->decay_momentum_b(_s, _sac);

            // Update floating masses in the box and evaluate
            _B.set_external_masses({_W, M_PION, M_JPSI, M_PION});
            _B.set_invariant_masses( _sac, _sab );
            _vBc = _B.eval_vector();
            
            // Dot products of momenta and box vector
            _q_dot_pb = q(x)*p_b(x) + q(y)*p_b(y) + q(z)*p_b(z);

            // go around the box accumulating the common couplings
            _C  = gy / sqrt(2.) * sqrt(M_Y*M_D1*M_D); 
            // Skip the D1 -> D* pi coupling which we include above
            _C *= G_PSI        * sqrt(M_JPSI*M_DSTAR*M_DSTAR);
            _C *= G1_PION      * sqrt(M_DSTAR * M_D);
        };

        // Assemble the vector containing all the box functions
        // p_a = jpsi,  p_b = pi_2,  p_c = pi_1
        inline complex q(cartesian_index i)
        {
            return 2.*_mpb*(_vBc[2]*p_c(i) - _vBc[3]*p_b(i)) + _vBc[0]*(p_c(i) - p_b(i));
        };

        // The box vector needs to be contracted with the jpsi vertex
        // also multiply by common factors
        inline complex M(cartesian_index j, cartesian_index k)
        {
            return _C*(q(k)*p_b(j) - q(j)*p_b(k) - delta(j,k)*_q_dot_pb);
        };

        // Scalar prefactors 
        complex _C = 0.;

        // Momenta
        double _mpc = 0., _mpb = 0.;  

        // Dot products of momenta
        complex _q_dot_pb = 0., _q_dot_pc = 0;

        // Box functions
        // Fixed masses for all the particles in the box
        box _B;
        std::array<double, 4> _internal = {M_D1, M_D, M_DSTAR, M_DSTAR}; 
        std::array<double, 4> _external = {4.2, M_PION, M_JPSI, M_PION};

        // Array to store the vector components of the box integral
        std::array<complex,4> _vBc, _vBb;

        // We need access to the D1D molecular nature of the Y state
        molecule _Y;
    };

};

#endif