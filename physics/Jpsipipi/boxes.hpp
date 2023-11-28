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
#include "Y(4260).hpp"

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
                result += D1_coupling(i, j) * M(j, k);
                // particle b couples to the D1 vertex
                _pi1 = particle::b;
                result += D1_coupling(i, j) * M(j, k);
            };
            return  _C * result / sqrt(2.);
        };

        protected:

        // Save which particle is considered pi1, by default this is particle c
        particle _pi1 = particle::c;

        // D1 -> D* pi coupling
        inline complex D1_coupling(cartesian_index i, cartesian_index j)
        {
            double swave = H1_S*sqrt(modp1()*modp1() + M_PION*M_PION) * delta(i,j);
            double dwave = H1_D*(3.*p1(i)*p1(j) - delta(i,j)*modp1()*modp1());

            return sqrt(M_D10*M_DSTAR0)*(swave + dwave);
        };

        inline void recalculate()
        {
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
            _B.set_invariant_masses(permute(_schan), permute(_tchan));
            _vB[1] = _B.eval_vector();
            _q_dot_p[1] =  q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);
        
            // go around the box accumulating the common couplings
            _C  = _Y->coupling()/sqrt(2.) * sqrt(_Y->mass()*M_D1*M_D0); 
            // Skip the D1 -> D* pi coupling which we include above
            _C *= _gjpsi;
            _C *= _gpi;
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
        inline double modp1(){ return (_pi1 == c) ? _mpc : _mpb; };
        inline double modp2(){ return (_pi1 == c) ? _mpb : _mpc; };
        inline double p1(cartesian_index i){ return modp1() * phat_1(i); };
        inline double p2(cartesian_index i){ return modp2() * phat_2(i); };

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
        double _C = 0., _gjpsi = 0., _gpi = 0.; 

        std::array<double,3> _decay_masses, _loop_masses;

        // Box functions
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
            _schan = ab, _tchan = ac;

            _decay_masses = {M_PION,  M_JPSI,   M_PION};
            _loop_masses  = {M_D0,    M_DSTAR0, M_DSTAR0};

            _gjpsi =  2.*G_PSI * sqrt(M_JPSI*M_DSTAR0*M_DSTAR0);
            _gpi   =  G1_PION; 
        };

        protected:

        // Things that specify the box
        // p1 == pi_1, p2 == pi_2, p3 == jpsi
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();

            complex pi1_piece = + (vecB[0] + 2.*vecB[2]);
            complex pi2_piece = - (vecB[0] + 2.*vecB[3]);
            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return p2(j)*q(k) - p2(k)*q(j) - delta(j, k)*qdotp2();
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

            _decay_masses = {M_JPSI,  M_PION,   M_PION};
            _loop_masses  = {M_D0,    M_DSTAR0, M_DSTAR0};

            _gjpsi = 2.*G_PSI * sqrt(M_JPSI*M_DSTAR0*M_D0); 
            _gpi   = G2_PION; 
        };

        protected:

        // Things that specify the box
        // p1 == pi_1, p2 == pi_2, p3 == jpsi
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();
            complex pi1_piece = - (2.*vecB[3] + 2.*vecB[2] + vecB[0]);
            complex pi2_piece = - (2.*vecB[3]              + vecB[0]);
            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return p2(k)*q(j) - delta(j, k)*qdotp2();
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

            _decay_masses = {M_JPSI,  M_PION,   M_PION};
            _loop_masses  = {M_D0,    M_D0,     M_DSTAR0};

            _gjpsi = 2.*G_PSI * sqrt(M_JPSI*M_D0*M_D0); 
            _gpi   = G1_PION; 
        };

        protected:

        // Things that specify the box
        // p1 == pi_1, p2 == pi_2, p3 == jpsi
        inline complex q(cartesian_index i)
        {
            auto vecB = vB();
            complex pi1_piece = - (2.*vecB[3] + 2.*vecB[2] + vecB[0]);
            complex pi2_piece = - (2.*vecB[3]              + vecB[0]);
            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return p2(j)*q(k);
        };
    };
};

#endif