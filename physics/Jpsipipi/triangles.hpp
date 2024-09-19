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
        _T(triangle::kLoopTools), _TNR(triangle::kNonrelativistic)
        {
            // Set up Zc propagator
            _Zc->set_parameters({3.9, 50.E-3, 4.6865/sqrt(M_DSTAR0*M_D0*3.9)});
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index k)
        {
            // print("right peak");
            // _pi1 = particle::c;
            // for (auto i : C_INDICES)
            // {
            //     for (auto j : C_INDICES)
            //     {
            //         print(i+1, j+1, M(i, j));
            //     }
            // }
            // line();

            // print("left peak");
            // _pi1 = particle::b; 
            // for (auto i : C_INDICES)
            // {
            //     for (auto j : C_INDICES)
            //     {
            //         print(i+1, j+1, M(i, j));
            //     }
            // }
            // exit(1);


            // _pi1 = particle::c; // Right peak
            // return M(i, k);

            complex result = 0.;
            for (auto j : C_INDICES)
            {
                // // particle c couples to the D1 vertex
                // _pi1 = particle::c; // Right peak
                // result += D1_coupling(i, j) * M(j, k);

                // particle b couples to the D1 vertex
                _pi1 = particle::b; // Left peak
                result += D1_coupling(i, j) * M(j, k);
            };

            return  _C * result / sqrt(2.);
        };

        protected:

        bool _printed = false;
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
            // ------------------------------------------------------
            // Internal masses of triangle 1 are common to both
            // pion configurations
            _TNR.set_internal_masses({M_DSTAR, M_D1, M_D}); 
            _TNR.add_width(2, W_D1); 
            _TNR.set_external_masses({_W, sqrt(_sab), M_PION});
            _T1[0] = _TNR.eval() * _Zc->propagator(_sab);
            _TNR.set_external_masses({_W, sqrt(_sac), M_PION});
            _T1[1] = _TNR.eval() * _Zc->propagator(_sac);

            // ------------------------------------------------------
            // Set up triangle 2
            _pi1 = particle::c; // c couples to D1
            _T.set_internal_masses({_loop_masses[0], _loop_masses[1], _loop_masses[2]});
            _T.set_external_masses({sqrt(_sab), M_JPSI, M_PION});
            _vT2[0]     = _T.eval_vector();
            _q_dot_p[0] = q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);

            // ------------------------------------------------------
            // Swap pions 
            _pi1 = particle::b; // b couples to D1
            _T.set_external_masses({sqrt(_sac), M_JPSI, M_PION});
            _vT2[1]     = _T.eval_vector();
            _q_dot_p[1] = q(x)*p2(x) + q(y)*p2(y) + q(z)*p2(z);

            // ------------------------------------------------------
            // Couplings at the vertices of the first triangle
            _C  = _Y->coupling()/sqrt(2.) * sqrt(_Y->mass()*M_D1*M_D0);     // Y -> D1 D
            // Skip the D1 -> D* pi coupling which we include above
            _C *= _Zc->coupling() * sqrt(M_D0*M_DSTAR0*_Zc->mass()); // D* D -> Z 
            // Gather couplings from the vertices of the second triangle
            _C *= _Zc->coupling() * sqrt(M_D0*M_DSTAR0*_Zc->mass()); // Z decay vertex
            _C *= _gpi;                                              // Coupling at pion vertex
            _C *= _gjpsi;                                            // Coupling at psi vertex 
        };

        // Aliases for the pion momenta specifying which is pi1 and pi2
        inline double modp1(){ return (_pi1 == c) ? _mpc : _mpb; };
        inline double modp2(){ return (_pi1 == c) ? _mpb : _mpc; };
        inline double p1(cartesian_index i){ return modp1() * phat_1(i); };
        inline double p2(cartesian_index i){ return modp2() * phat_2(i); };

        // Returns T1 for current permutation of pions
        inline complex T1(){ return (_pi1 == particle::c) ? _T1[0] : _T1[1]; };

        // Returns vector decomposition of T2 for current permutation of pions
        inline std::array<complex,3> vT2(){ return (_pi1 == particle::c) ? _vT2[0] : _vT2[1]; };

        // Unlike the boxes, the triangle doesnt involve reordering the external masses so 
        // we have a common q
        // p1 == pi_1, p2 == pi_2, p3 == jpsi
        inline complex q(cartesian_index i)
        {
            auto vecT2 = vT2();
            complex pi1_piece = - (2.* vecT2[1] + 2.*vecT2[2] + vecT2[0]);
            complex pi2_piece = - (2.* vecT2[2]               + vecT2[0]);
            return pi1_piece * p1(i) + pi2_piece * p2(i);
        };
        inline complex qdotp2(){ return (_pi1 == particle::c) ? _q_dot_p[0] : _q_dot_p[1]; };
        // Assemble the vector containing all the triangle functions functions
        // This will be contracted with the jpsi vertex
        // also multiply by common factors
        virtual inline complex M(cartesian_index j, cartesian_index k) = 0;

        // Scalar prefactors 
        double _C     = 0.;

        // Couplings to be specified by the specific implementation
        double _gjpsi = 0., _gpi = 0.; 

        std::array<double,3> _loop_masses;

        // Triangle functions, this is used to evaluate both triangles
        triangle _T, _TNR;

        // Cached quantities related to the triangle
        // index 0 couples particle c to the D-wave vertex while 1 is particle b
        std::array<complex,2>     _q_dot_p, _T1;  // Scalar Y -> Zc pi triangle
        std::array<std::array<complex,3>,2> _vT2; // Vector Zc -> Jpsi pi triangle

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
            _gjpsi       = 2.*G_PSI*sqrt(M_JPSI * M_DSTAR0 * M_D0); 
            _gpi         = G2_PION;
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
            _gjpsi       = 2.*G_PSI * sqrt(M_JPSI  * M_D0 * M_D0); 
            _gpi         = G1_PION;
        };

        protected:

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return /* T1() * */ p2(j)*q(k);
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
            _gjpsi       = 2.*G_PSI * sqrt(M_JPSI  * M_DSTAR0 * M_DSTAR0); 
            _gpi         = G1_PION;
        };

        protected:

        inline complex M(cartesian_index j, cartesian_index k)
        {
            return /* T1() * */ (p2(j)*q(k) - p2(k)*q(j) - delta(j, k)*qdotp2());
        };
    };    
};

#endif