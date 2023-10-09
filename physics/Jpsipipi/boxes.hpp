// Specific implementations of box amplitudes relevant for the Jpsi pi pi final state

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
                // particle c couples to the D1 vertex
                _pi1 = particle::c;
                result += d_wave(i, j) * M(j, k);

                // particle b couples to the D1 vertex
                _pi1 = particle::b;
                result += d_wave(i, j) * M(j, k);
            };
            return  _C * result / sqrt(2.);
        };

        protected:

        // Save which particle is considered pi1, by default this is particle c
        particle _pi1 = particle::c;

        // D1 -> D* pi coupling
        inline complex d_wave(cartesian_index i, cartesian_index j)
        {
            return sqrt(M_D1*M_DSTAR) * (H1_D*(3.*p1(i)*p1(j) - delta(i,j)*modp()*modp()) + H1_S*delta(i,j));
        };

        inline void recalculate()
        {
            // Y mass and couplings
            double gy   = _Y->molecular_coupling();
            double M_Y = sqrt(_s);

            // Update floating masses in the box and evaluate
            _B.set_external_masses({_W,   _decay_masses[0], _decay_masses[1], _decay_masses[2]});
            _B.set_internal_masses({M_D1, _loop_masses[0],  _loop_masses[1],  _loop_masses[2] });
            _B.add_width(0, W_D1);

            // Box vector related things depend on which particle is particle 1
            _pi1 = particle::c;
            _B.set_invariant_masses(translate(_schan), translate(_tchan));
            _vB[0] = _B.eval_vector();
            _q_dot_p[0] =  q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);
            
            _pi1 = particle::b;
            _B.set_invariant_masses(permute(_schan),   permute(_tchan));
            _vB[1] = _B.eval_vector();
            _q_dot_p[1] =  q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);
        
            // go around the box accumulating the common couplings
            _C  = gy / sqrt(2.) * sqrt(M_Y*M_D1*M_D); 
            // Skip the D1 -> D* pi coupling which we include above
            _C *= _gjpsi;
            _C *= G1_PION       * sqrt(M_DSTAR * M_D);
        };

        // Methods related to handling invariant masses for the box
        // translate gives the specified subchannel invariant mass
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
        // Permute takes the given subchannel swaps paritlces b <-> c 
        // and returns the invariant mass of the new subchannel
        inline double permute(subchannel sig)
        {
            switch (sig)
            {
                case ab: return _sac;
                case bc: return _sbc;
                case ac: return _sab;
            };
        };

        // Aliases for the pion momenta specifying which is pi1 and pi2
        inline double modp(){ return (_pi1 == c) ? _mpc : _mpb; };
        inline complex p1(cartesian_index i){ return modp() * phat_1(i); };
        inline complex p2(cartesian_index i){ return modp() * phat_2(i); };

        // Vector decomposition of the Box
        inline std::array<complex,4> vB(){ return (_pi1 == particle::c) ? _vB[0] : _vB[1]; };

        // Assemble the vector containing all the box functions
        // p_a = jpsi,  p_b = pi_2,  p_c = pi_1
        virtual inline complex q(cartesian_index i) = 0;
        inline complex qdotp2(){ return (_pi1 == particle::c) ? _q_dot_p[0] : _q_dot_p[1]; };

        // The box vector needs to be contracted with the jpsi vertex
        // also multiply by common factors
        virtual inline complex M(cartesian_index j, cartesian_index k) = 0;

        // Scalar prefactors 
        double _C     = 0.;
        double _gjpsi = 0.; 

        std::array<double,3> _decay_masses;
        std::array<double,3> _loop_masses;

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
    // Specific boxes follow notation in Leon's note

    // Box with J/psi D* D* coupling
    class box_I : public generic_box
    {
        public: 

        box_I(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box I")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = ac, _tchan = ab;

            _decay_masses = {M_PION, M_JPSI,  M_PION};
            _loop_masses  = {M_D,    M_DSTAR, M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_DSTAR*M_DSTAR); 
        };

        protected:

        // Things that specify the box
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();
            complex pi1_piece =  2.*vecB[2] + vecB[0];
            complex pi2_piece = -2.*vecB[3] - vecB[0];

            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return q(k)*p2(j) - q(j)*p2(k) - delta(j,k)*qdotp2();
        };
    };

    // Box with J/psi D D* coupling
    class box_II : public generic_box
    {
        public: 

        box_II(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box II")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = bc, _tchan = ab;

            _decay_masses = {M_JPSI, M_PION,  M_PION};
            _loop_masses  = {M_D,    M_DSTAR, M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_DSTAR*M_D); 
        };

        protected:

        // Things that specify the box
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();
            complex pi1_piece =  2.*(vecB[3] + vecB[2] + vecB[0]);
            complex pi2_piece =  2.* vecB[3] + vecB[0];

            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return q(j)*p2(k) - delta(j,k)*qdotp2();
        };
    };

    // Box with J/psi D D coupling
    class box_III : public generic_box
    {
        public: 

        box_III(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box III")
        : generic_box(key, xkinem, Y, id)
        {
            _schan = bc, _tchan = ab;

            _decay_masses = {M_JPSI, M_PION,  M_PION};
            _loop_masses  = {M_D,    M_D,     M_DSTAR};

            _gjpsi = G_PSI * sqrt(M_JPSI*M_D*M_D); 
        };

        protected:

        // Things that specify the box
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();
            complex pi1_piece =  2.*(vecB[3] + vecB[2] + vecB[0]);
            complex pi2_piece =  2.* vecB[3] + vecB[0];

            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return q(k)*p2(j);
        };
    };
};

#endif