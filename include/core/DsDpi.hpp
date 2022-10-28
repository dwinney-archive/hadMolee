// Specific implementations of amplitudes relevant for the D* D pi final state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DSDPI
#define DSDPI

#include "reaction_kinematics.hpp"
#include "breit_wigner.hpp"
#include "amplitude.hpp"
#include "constants.hpp"
#include "triangle.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Contact-like interaction in S-wave
    // This contains the exotic Zc channel, we multiply by a linear polynomial and fit the coefficients

    class DsDpi_swave : public amplitude
    {
        // -----------------------------------------------------------------------
        public:
        
        DsDpi_swave(reaction_kinematics * xkinem, D1D_molecule * Y, std::string id = "DsDpi_swave")
        : amplitude(xkinem, Y, 2, "DsDpi_swave", id)
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline std::complex<double> reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // IF our kinematics has been changed, recalculate relevant quantities
            if (updated()) recalculate();

            // Being S-wave the s-wave strength gets multiplied by a delta-function
            return _AS * s_wave(i,j);
        };

        // Set the two parameters
        // These are the coefficients of the linear polynomial multiplying the s-wave strength
        inline void set_parameters(std::vector<double> par)
        {
            check_nParams(par);

            _a = par[0];
            _b = par[1];
        };

        // -----------------------------------------------------------------------
        private:

        // Parameterize with two real polynomial coefficients
        double _a = A_QQ2016, _b = B_QQ2016;

        // S-wave coupling strength
        std::complex<double> _AS; 

        // This channel recieves contribution from the Z(3900)
        DsD_molecule _Zc;

        // S-wave amplitude is easy, because its only the Z propagator that needs to be calcualted
        inline void recalculate()
        {
            // sab is assumed to be DsD channel
            _AS = _a * (_sab + _b) * _Zc.propagator(_sab);
        };
    };

    // ---------------------------------------------------------------------------
    // Tree level transition through intermediate D1
    // No free parameters since this assume to be known background

    class DsDpi_tree : public amplitude
    {
        // -----------------------------------------------------------------------
        public:
        
        DsDpi_tree(reaction_kinematics * xkinem, D1D_molecule * Y, std::string id = "DsDpi_tree")
        : amplitude(xkinem, Y, 0, "DsDpi_tree", id),  _Y(Y), _D1(M_D1, W_D1)
        {};

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline std::complex<double> reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            if (updated()) recalculate();

            return _AS * s_wave(i,j) + _AD * d_wave(i,j);
        };

        // -----------------------------------------------------------------------
        private:

        inline void recalculate()
        {
            double p_pion  = _kinematics->decay_momentum_c(_s, _sab);
            double p2_pion = p_pion*p_pion;

            // sab is assumed to be DsD channel
            _AS = XI * _hS * sqrt(p2_pion + _mc2);
            _AD = XI * _hD * p2_pion;

            // Multiply by the propagator of the D1 and y coupling
            _AS *= (_Y->molecular_coupling() / sqrt(2.)) * _D1.eval(_sac);
            _AD *= (_Y->molecular_coupling() / sqrt(2.)) * _D1.eval(_sac);
        };

        double _hS = 0., _hD = HPRIME_UPPER / F_PION;    // Coupling constants to S and D wave interactions
        std::complex<double> _AS, _AD;               // Energy dependent S and D wave strengths

        // This channel has a D1 resonance in the Ds pi channel
        relativistic_BW _D1;

        // And we explicitly require information of molecular nature of Y 
        // so we save a D1D_molecule version of the pointer, not just charmoniumlike
        D1D_molecule *_Y;
    };

    // ---------------------------------------------------------------------------
    // Transition involving D1DsD triangle and Z meson
    // No free parameters since this assume to be known background

    class DsDpi_triangle : public amplitude
    {
        // -----------------------------------------------------------------------
        public:
        
        // Here we can choose whether we want a nonrelativistic triangle or the relativistic version
        // We default to the nonrel version
        DsDpi_triangle(reaction_kinematics * xkinem, D1D_molecule * Y, std::string id = "DsDpi_triangle")
        : amplitude(xkinem, Y, 0, "DsDpi_triangle", id), _Y(Y), _T(new nonrelativistic_triangle(_external, _internal))
        {};

        // Destructor needs to clean up the new pointer we made
        ~DsDpi_triangle()
        {
            delete _T;
        };

        // The reduced amplitude corresponds to the S-wave contact-like interaction
        // however it receives contributions from the propagator of the Z meson
        inline std::complex<double> reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            // Check if we need to recalculate strengths
            if (updated()) recalculate();

            // Being D-wave we get the appropriate projector
            // Assume pion (particle c) defined the +z direction
            return _AS * s_wave(i,j) + _AD * d_wave(i,j);
        };

        // Switch the evaluation of the triangle to the relativistic version
        inline void use_relativistic(bool x = true)
        {
            // Delete the existing pointer
            delete _T;

            // and reallocate
            if (x) { _T = new relativistic_triangle   (_external, _internal); } 
            else   { _T = new nonrelativistic_triangle(_external, _internal); };
        };

        // -----------------------------------------------------------------------
        private:
        
        inline void recalculate()
        {
            double p_pion = _kinematics->decay_momentum_c(_s, _sab);

            // sab is assumed to be DsD channel
            _AS = XI * _hS * sqrt(p_pion*p_pion + _mc2);
            _AD = XI * _hD * p_pion*p_pion;

            // Update arguments with floating Y and Z meson masses
            _T->set_external_masses({_W, sqrt(_sab), M_PION}); 
            
            // Multiply by the propagator of the Z and triangle function
            double z = _Zc.molecular_coupling();
            
            // The normalization is matched to Qiang's paper
            complex<double> T = _internal[0] * _internal[1] * _internal[2] * _T->eval();

            _AS *= - (_Y->molecular_coupling() / sqrt(2.)) * z*z * T * _Zc.propagator(_sab);
            _AD *= - (_Y->molecular_coupling() / sqrt(2.)) * z*z * T * _Zc.propagator(_sab);
        };

        double _hS = 0., _hD = HPRIME_UPPER / F_PION;  // D1 -> D*pi coupling for the S-wave and the D-wave
        std::complex<double> _AS, _AD;                 // Energy dependent S and D wave strengths

        // On-shell masses involved in the triangle
        std::array<double,3> _internal = {M_DSTAR, M_D1,     M_D};
        std::array<double,3> _external = {M_Y4260, M_ZC3900, M_PION};

        // Default to using nonrelativistic version of the triangle
        triangle *_T;

        // Z meson resonance in this diagram
        DsD_molecule _Zc;

        // And we explicitly require information of molecular nature of Y 
        // so we save a D1D_molecule version of the pointer, not just charmoniumlike
        D1D_molecule *_Y;
    };
};

#endif