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
    // Generic box amplitude which is automatically symmetrices with regards to the two pions
    // Just need to specify the channels and parcitles propagating in the loop 
    class generic_box : public amplitude_base
    {
        public:
        
        generic_box(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "generic_box")
        : amplitude_base(key, xkinem, Y, 0, "generic_box", id), 
        _B(box::kLoopTools), _Y(get_molecular_component(Y))
        {};

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index k)
        {
            complex result = 0.;
            for (auto j : C_INDICES)
            {
                // pi1 couples to the D1 vertex
                result += (H1_D*(3.*p(c, i)*p(c, j) - delta(i,j)*_mpc*_mpc) + H1_S*delta(i,j)) * M(b, j, k);

                // pi2 couples to the D1 vertex
                result += (H1_D*(3.*p(b, i)*p(b, j) - delta(i,j)*_mpb*_mpb) + H1_S*delta(i,j)) * M(c, j, k);
            };
            return  _C * result / sqrt(2.);
        };

        protected:

        inline void recalculate()
        {
            // Y mass and couplings
            double gy   = _Y->molecular_coupling();
            double M_Y = sqrt(_s);

            // Update floating masses in the box and evaluate
            _B.set_external_masses({_W,   _decay_masses[0],    _decay_masses[1],    _decay_masses[2]});
            _B.set_internal_masses({M_D1, _internal_masses[0], _internal_masses[1], _internal_masses[2]});
            _B.add_width(0, W_D1);

            _B.set_invariant_masses(translate(_schan), translate(_tchan));
            _vB[0] = _B.eval_vector();
            
            _B.set_invariant_masses(permute(_schan),   permute(_tchan));
            _vB[1] = _B.eval_vector();
        
            // Dot products of momenta and box vector
            _q_dot_p[0] =  q(b, x)*p_b(x) + q(b, y)*p_b(y) + q(b, z)*p_b(z);
            _q_dot_p[1] =  q(c, x)*p_c(x) + q(c, y)*p_c(y) + q(c, z)*p_c(z);

            // go around the box accumulating the common couplings
            _C  = gy / sqrt(2.) * sqrt(M_Y*M_D1*M_D); 
            // Skip the D1 -> D* pi coupling which we include above
            _C *= _gjpsi;
            _C *= G1_PION       * sqrt(M_DSTAR * M_D);
        };

        //;
        subchannel _schan, _tchan;
        inline double translate(subchannel sig)
        {
            switch (sig)
            {
                case ab: return _sab;
                case bc: return _sbc;
                case ac: return _sac;
            };
        };
        inline double permute(subchannel sig)
        {
            switch (sig)
            {
                case ab: return _sac;
                case bc: return _sbc;
                case ac: return _sab;
            };
        };

        inline particle other_pion(particle pi1)
        {
            return (pi1==b) ? c : b;
        };

        inline std::array<complex,4> vB(particle pi2)
        {
            int n = (pi2 == b) ? 0 : 1;
            return _vB[n];
        };

        inline complex qdotp(particle pi2)
        {
            return (pi2 == b) ? _q_dot_p[0] : _q_dot_p[1];
        };

        // Assemble the vector containing all the box functions
        // p_a = jpsi,  p_b = pi_2,  p_c = pi_1
        virtual inline complex q(particle pi2, cartesian_index i) = 0;

        // The box vector needs to be contracted with the jpsi vertex
        // also multiply by common factors
        virtual inline complex M(particle pi2, cartesian_index j, cartesian_index k) = 0;

        // Scalar prefactors 
        double _C     = 0.;
        double _gjpsi = 0.; 

        std::array<double,3> _decay_masses;
        std::array<double,3> _internal_masses;

        // Box functions
        // Fixed masses for all the particles in the box
        box _B;

        // Array to store the vector components of the box integral
        std::array<complex,2> _q_dot_p;
        std::array<std::array<complex,4>,2> _vB;

        // We need access to the D1D molecular nature of the Y state
        molecule _Y;
    };

    // ---------------------------------------------------------------------------
    // Specific boxes follows notation in Leon's note
    
    class Box_I : public generic_box
    {
        public: 

        Box_I(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box I")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = ac, _tchan = ab;

            _decay_masses    = {M_PION, M_JPSI,  M_PION};
            _internal_masses = {M_D,    M_DSTAR, M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_DSTAR*M_DSTAR); 
        };

        protected:

        // Things that specify the box
        inline complex q(particle pi2, cartesian_index i)
        {
            auto vecB = vB(pi2);
            complex pi1_piece =  2.*vecB[2] + vecB[0];
            complex pi2_piece = -2.*vecB[3] - vecB[0];

            return pi1_piece * p(other_pion(pi2), i) + pi2_piece * p(pi2, i);
        };

        inline complex M(particle pi2, cartesian_index j, cartesian_index k)
        {
            return q(pi2, k)*p(pi2, j) - q(pi2, j)*p(pi2, k) - delta(j,k)*qdotp(pi2);
        };
    };


    class Box_II : public generic_box
    {
        public: 

        Box_II(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box II")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = bc, _tchan = ab;

            _decay_masses    = {M_JPSI, M_PION,  M_PION};
            _internal_masses = {M_D,    M_DSTAR, M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_DSTAR*M_D); 
        };

        protected:

        // Things that specify the box
        inline complex q(particle pi2, cartesian_index i)
        {
            auto vecB = vB(pi2);
            complex pi1_piece =  2.*(vecB[3] + vecB[2] + vecB[0]);
            complex pi2_piece =  2.* vecB[3] + vecB[0];

            return pi1_piece * p(other_pion(pi2), i) + pi2_piece * p(pi2, i);
        };

        inline complex M(particle pi2, cartesian_index j, cartesian_index k)
        {
            return q(pi2, j)*p(pi2, k) - delta(j,k)*qdotp(pi2);
        };
    };

    class Box_III : public generic_box
    {
        public: 

        Box_III(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box III")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = bc, _tchan = ab;

            _decay_masses    = {M_JPSI, M_PION,  M_PION};
            _internal_masses = {M_D,    M_D,     M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_DSTAR*M_D); 
        };

        protected:

        // Things that specify the box
        inline complex q(particle pi2, cartesian_index i)
        {
            auto vecB = vB(pi2);
            complex pi1_piece =  2.*(vecB[3] + vecB[2] + vecB[0]);
            complex pi2_piece =  2.* vecB[3] + vecB[0];

            return pi1_piece * p(other_pion(pi2), i) + pi2_piece * p(pi2, i);
        };

        inline complex M(particle pi2, cartesian_index j, cartesian_index k)
        {
            return q(pi2, k)*p(pi2, j);
        };
    };
};

#endif