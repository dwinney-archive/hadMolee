// Specific implementations of amplitudes relevant for the Jpsi pi pi final state

// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef JPSIPIPI_TRIANGLES_HPP
#define JPSIPIPI_TRIANGLES_HPP

#include <memory>

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "triangle.hpp"
#include "Y(4260).hpp"

namespace hadMolee::Jpsipipi
{
    // Generic triangle amplitude which is automatically symmetrices with regards to the two pions
    // Just need to specify the channels and parcitles propagating in the loop 
    class generic_triangle : public amplitude_base
    {
        public:
        
        generic_triangle(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "generic_box")
        : amplitude_base(key, xkinem, Y, 0, "generic_triangle", id), 
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
            return sqrt(M_D1*M_DSTAR)*(H1_D*(3.*p1(i)*p1(j) - delta(i,j)*modp()*modp()) + H1_S*delta(i,j));
        };

        inline void recalculate()
        {
            // Y mass and couplings
            double gy  = _Y->coupling();
            double M_Y = sqrt(_s);
            double gz  = _Zc->coupling();
            double M_Z = M_ZC3900;

            // Update floating masses in triangle 1 and evaluate
            // Also include the Zc propagator here
            _T.set_external_masses({_W, sqrt(_sab), M_PION});
            _T.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _T.add_width(2, W_D1); 
            _T1[0] = _T.eval() * _Zc->propagator(_sab);
            // Swap pions 
            _T.set_external_masses({_W, sqrt(_sac), M_PION});
            _T1[1] = _T.eval() * _Zc->propagator(_sac);

            // Couplings at the vertices of the first triangle
            _C  = gy/sqrt(2.) * sqrt(M_Y*M_D1*M_D);    // Y -> D1 D
            // Skip the D1 -> D* pi coupling which we include above
            _C *= gz          * sqrt(M_D*M_DSTAR*M_Z); // D* D -> Z 

            // Set up triangle 2
            // This will depend on the running masses in specific triangle diagram
            _pi1 = particle::c;
            _T.set_external_masses({sqrt(_sab), M_JPSI, M_PION});
            _T.set_internal_masses({_loop_masses[0], _loop_masses[1], _loop_masses[2]});
            _vT2[0]     = _T.eval_vector();
            _q_dot_p[0] = q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);
            // Swap pions 
            _pi1 = particle::b;
            _T.set_external_masses({sqrt(_sac), M_JPSI, M_PION});
            _vT2[1]     = _T.eval_vector();
            _q_dot_p[1] = q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);

            // Gather couplings from the vertices of the second triangle
            _C *= gz     * sqrt(M_D*M_DSTAR*M_Z); // Z decay vertex
            _C *= _gpi;                           // Coupling at pion vertex
            _C *= _gjpsi;                         // Coupling at psi vertex 
        };

        // Aliases for the pion momenta specifying which is pi1 and pi2
        inline double modp(){ return (_pi1 == c) ? _mpc : _mpb; };
        inline complex p1(cartesian_index i){ return modp() * phat_1(i); };
        inline complex p2(cartesian_index i){ return modp() * phat_2(i); };

        // Returns T1 for current permutation of pions
        inline complex T1(){ return (_pi1 == particle::c) ? _T1[0] : _T1[1]; };

        // Returns vector decomposition of T2 for current permutation of pions
        inline std::array<complex,3> vT2(){ return (_pi1 == particle::c) ? _vT2[0] : _vT2[1]; };

        // Assemble the vector containing all the box functions
        inline complex qdotp2(){ return (_pi1 == particle::c) ? _q_dot_p[0] : _q_dot_p[1]; };
        
        // Unlike the boxes, the triangle doesnt involve reordering the external masses so 
        // we have a common q
        // p1 == pi_1, p2 == pi_2, p3 == jpsi
        inline complex q(cartesian_index i)
        {
            auto vecT2 = vT2();

            // 2 l - (p2 + p1)
            // -> 2 (C1 p2 + C2 (p2 + p3)) - C0 (p1 + p2)
            // =  2 (C1 p2 - C2 p1) - C0 (p1 + p2)
            complex pi1_piece = - 2.* vecT2[1] + vecT2[0];
            complex pi2_piece = + 2.* vecT2[2] - vecT2[0];
            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };
        // The box vector needs to be contracted witlsh the jpsi vertex
        // also multiply by common factors
        virtual inline complex M(cartesian_index j, cartesian_index k) = 0;

        // Scalar prefactors 
        double _C     = 0.;

        // Couplings to be specified by the specific implementation
        double _gjpsi = 0., _gpi = 0.; 

        std::array<double,3> _loop_masses;

        // Triangle functions, this is used to evaluate both triangles
        triangle _T;

        // Cached quantities related to the triangle
        // index 0 couples particle c to the D-wave vertex while 1 is particle b
        std::array<complex,2>               _T1;  // Scalar Y -> Zc pi triangle
        std::array<std::array<complex,3>,2> _vT2; // Vector Zc -> Jpsi pi triangle
        std::array<complex,2> _q_dot_p;

        // We need access the molecular nature of the Y and Zc state
        molecule _Y, _Zc;
    };

    // ---------------------------------------------------------------------------
    // Specific triangles follow notation in Leon's note

    // Triangle with Jpsi D D* vertex
    class triangle_I : public generic_triangle
    {
        public:

        triangle_I(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Triangle I")
        : generic_triangle(key, xkinem, Y, id)
        {
            _loop_masses = {M_D, M_DSTAR, M_DSTAR};
            _gjpsi       = G_PSI   * sqrt(M_JPSI  * M_DSTAR * M_D); 
            _gpi         = G2_PION * sqrt(M_DSTAR * M_DSTAR);
        };

        protected:

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return T1() * (p2(k)*q(j) - delta(j, k)*qdotp2());
        };
    };  

    // Triangle with Jpsi D D vertex
    class triangle_II : public generic_triangle
    {
        public:

        triangle_II(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Triangle II")
        : generic_triangle(key, xkinem, Y, id)
        {
            _loop_masses = {M_D, M_DSTAR, M_D};
            _gjpsi       = G_PSI   * sqrt(M_JPSI  * M_D * M_D); 
            _gpi         = G1_PION * sqrt(M_DSTAR * M_D);
        };

        protected:

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return T1() * p2(j)*q(k);
        };
    };  

    // Triangle with Jpsi D* D* vertex
    class triangle_III : public generic_triangle
    {
        public:

        triangle_III(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Triangle III")
        : generic_triangle(key, xkinem, Y, id)
        {
            _loop_masses = {M_DSTAR, M_D, M_DSTAR};
            _gjpsi       = G_PSI   * sqrt(M_JPSI  * M_DSTAR * M_DSTAR); 
            _gpi         = G1_PION * sqrt(M_DSTAR * M_D);
        };

        protected:

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return T1() * (p2(j)*q(k) - p2(k)*q(j) - delta(j,k)*qdotp2());
        };
    };    
};

#endif