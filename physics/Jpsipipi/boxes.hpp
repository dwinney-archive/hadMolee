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
        {};

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index k)
        {
            complex result = 0.;
            for (auto j : C_INDICES)
            {
                // pi1 couples to the D1 vertex
                result += (H1_D*(3.*p_c(i)*p_c(j) - delta(i,j)*_mpc*_mpc) + H1_S*delta(i,j)) * M_b(j, k);

                // pi2 couples to the D1 vertex
                result += (H1_D*(3.*p_b(i)*p_b(j) - delta(i,j)*_mpb*_mpb) + H1_S*delta(i,j)) * M_c(j, k);
            };
            return  result / sqrt(2.);
        };

        private:

        inline void recalculate()
        {
            // Y mass and couplings
            double gy   = _Y->molecular_coupling();
            double M_Y = sqrt(_s);

            // Update floating masses in the box and evaluate
            _B.set_external_masses({_W, M_PION, M_JPSI, M_PION});
            _B.set_internal_masses({M_D1, M_D, M_DSTAR, M_DSTAR});
            _B.add_width(0, W_D1);

            _B.set_invariant_masses( _sac, _sab );
            _vB[0] = _B.eval_vector();
            
            _B.set_invariant_masses( _sab, _sac );
            _vB[1] = _B.eval_vector();
        
            // Dot products of momenta and box vector
            _q_dot_p[0] =  q_b(x)*p_b(x) + q_b(y)*p_b(y) + q_b(z)*p_b(z);
            _q_dot_p[1] =  q_c(x)*p_c(x) + q_c(y)*p_c(y) + q_c(z)*p_c(z);

            // go around the box accumulating the common couplings
            _C  = gy / sqrt(2.) * sqrt(M_Y*M_D1*M_D); 
            // Skip the D1 -> D* pi coupling which we include above
            _C *= G_PSI         * sqrt(M_JPSI*M_DSTAR*M_DSTAR);
            _C *= G1_PION       * sqrt(M_DSTAR * M_D);
        };

        // Assemble the vector containing all the box functions
        // p_a = jpsi,  p_b = pi_2,  p_c = pi_1
        inline complex q_b(cartesian_index i)
        {
            return 2.*(_vB[0][2]*p_c(i) - _vB[0][3]*p_b(i)) + _vB[0][0]*(p_c(i) - p_b(i));
        };

        inline complex q_c(cartesian_index i)
        {
            return 2.*(_vB[1][2]*p_b(i) - _vB[1][3]*p_c(i)) + _vB[1][0]*(p_b(i) - p_c(i));
        };

        // The box vector needs to be contracted with the jpsi vertex
        // also multiply by common factors
        inline complex M_b(cartesian_index j, cartesian_index k)
        {
            return _C*(q_b(k)*p_b(j) - q_b(j)*p_b(k) - delta(j,k)*_q_dot_p[0]);
        };
        inline complex M_c(cartesian_index j, cartesian_index k)
        {
            return _C*(q_c(k)*p_c(j) - q_c(j)*p_c(k) - delta(j,k)*_q_dot_p[1]);
        };

        // Scalar prefactors 
        complex _C = 0.;

        // Box functions
        // Fixed masses for all the particles in the box
        box _B;

        // Array to store the vector components of the box integral
        std::array<complex,2> _q_dot_p;
        std::array<std::array<complex,4>,2> _vB;

        // We need access to the D1D molecular nature of the Y state
        molecule _Y;
    };

};

#endif