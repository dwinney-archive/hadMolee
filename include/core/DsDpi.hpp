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

// ---------------------------------------------------------------------------
// Contact-like interaction in S-wave

class DsDpi_swave : public amplitude
{
    // -----------------------------------------------------------------------
    public:
    
    DsDpi_swave(reaction_kinematics * xkinem, D1D_molecule * Y, string id = "DsDpi_swave")
    : amplitude(xkinem, Y, 2, "DsDpi_swave", id)
    {};

    // The reduced amplitude corresponds to the S-wave contact-like interaction
    // however it receives contributions from the propagator of the Z meson
    inline complex<double> reduced_amplitude(int i, int j)
    {

        if (updated()) recalculate();

        // Being S-wave the s-wave strength gets multiplied by a delta-function
        return _AS * s_wave(i,j);
    };

    // Set the two parameters
    inline void set_parameters(vector<double> par)
    {
        check_nParams(par);

        _a = par[0];
        _b = par[1];
    };

    // -----------------------------------------------------------------------
    private:

    // Parameterize with two real polynomial coefficients
    double _a = S_A, _b = S_B;

    // S-wave coupling strength
    complex<double> _AS; 

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
    
    DsDpi_tree(reaction_kinematics * xkinem, D1D_molecule * Y, string id = "DsDpi_tree")
    : amplitude(xkinem, Y, 0, "DsDpi_tree", id),  _Y(Y), _D1(M_D1, W_D1)
    {};

    // The reduced amplitude corresponds to the S-wave contact-like interaction
    // however it receives contributions from the propagator of the Z meson
    inline complex<double> reduced_amplitude(int i, int j)
    {
        if (updated()) recalculate();

        // Being D-wave we get the appropriate projector
        // Assume pion (particle c) defined the +z direction
        return _AS * s_wave(i,j) + _AD * d_wave(i,j);
    };

    // -----------------------------------------------------------------------
    private:

    inline void recalculate()
    {
        double p_pion  = _kinematics->decay_momentum_c(_s, _sab);
        double p2_pion = p_pion*p_pion;

        // Can use the debug flag to turn on and off the S-wave coupling
        if (_debug == 1) {_hS = 0.; _hD = HP * sqrt(6.);}

        // sab is assumed to be DsD channel
        _AS = XI * (_hS / F_PION) * (1./sqrt(6.)) * sqrt(p2_pion + _mc2);
        _AD = XI * (_hD / F_PION) * (1./sqrt(6.)) * p2_pion;

        // Multiply by the propagator of the D1 and y coupling
        _AS *= _Y->molecular_coupling() * _D1.eval(_sac);
        _AD *= _Y->molecular_coupling() * _D1.eval(_sac);
    };

    double _hS = HP_S, _hD = HP_D;  // We have DsDpi coupling for the S-wave and the D-wave
    complex<double> _AS, _AD;       // Energy dependent S and D wave strengths

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
    DsDpi_triangle(reaction_kinematics * xkinem, D1D_molecule * Y, string id = "DsDpi_triangle")
    : amplitude(xkinem, Y, 0, "DsDpi_triangle", id), _Y(Y), _T(new nonrelativistic_triangle(_external, _internal))
    {};

    // Destructor needs to clean up the new pointer we made
    ~DsDpi_triangle()
    {
        delete _T;
    };

    // The reduced amplitude corresponds to the S-wave contact-like interaction
    // however it receives contributions from the propagator of the Z meson
    inline complex<double> reduced_amplitude(int i, int j)
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

        // Can use the debug flag to turn on and off the S-wave coupling
        if (_debug == 1) {_hS = 0.; _hD = HP * sqrt(6.) ;}

        // sab is assumed to be DsD channel
        _AS = XI * (_hS / F_PION) * (1./sqrt(6.)) * sqrt(p_pion*p_pion + _mc2);
        _AD = XI * (_hD / F_PION) * (1./sqrt(6.)) * p_pion*p_pion;

        _T->set_external_masses({_W, sqrt(_sab), M_PION}); // Update arguments with floating Y and Z meson masses
        
        // Multiply by the propagator of the Z and triangle function
        double z = _Zc.molecular_coupling();

        _AS *= - _Y->molecular_coupling() * z*z * _T->eval() * _Zc.propagator(_sab);
        _AD *= - _Y->molecular_coupling() * z*z * _T->eval() * _Zc.propagator(_sab);
    };

    double _hS = HP_S, _hD = HP_D;  // D1 -> D*pi coupling for the S-wave and the D-wave
    complex<double> _AS, _AD;       // Energy dependent S and D wave strengths

    // On-shell masses involved in the triangle
    array<double,3> _internal = {M_DSTAR, M_D1,     M_D};
    array<double,3> _external = {M_Y4260, M_ZC3900, M_PION};

    // Default to using nonrelativistic version of the triangle
    triangle *_T;

    // Z meson resonance in this diagram
    DsD_molecule _Zc;

    // And we explicitly require information of molecular nature of Y 
    // so we save a D1D_molecule version of the pointer, not just charmoniumlike
    D1D_molecule *_Y;
};

#endif