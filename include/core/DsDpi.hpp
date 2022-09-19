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

class DsDpi_swave : public amplitude
{
    // -----------------------------------------------------------------------
    public:
    
    DsDpi_swave(reaction_kinematics * xkinem, hadronic_molecule * Y, string id = "DsDpi_swave")
    : amplitude(xkinem, Y, 2, "DsDpi_swave", id)
    {};

    // The reduced amplitude corresponds to the S-wave contact-like interaction
    // however it receives contributions from the propagator of the Z meson
    inline complex<double> reduced_amplitude(int i, int j)
    {
        // sab is assumed to be DsD channel
        complex<double> A_S = _a * (_sab + _b) * Z.propagator(_W);

        // Being S-wave the s-wave strength gets multiplied by a delta-function
        return A_S * (i == j);
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
    double _a, _b;

    // This channel recieves contribution from the Z(3900)
    DsD_molecule Z;
};

// Tree level transition through intermediate D1
// No free parameters since this assume to be known background
class DsDpi_tree : public amplitude
{
    // -----------------------------------------------------------------------
    public:
    
    DsDpi_tree(reaction_kinematics * xkinem, hadronic_molecule * Y, string id = "DsDpi_tree")
    : amplitude(xkinem, Y, 0, "DsDpi_tree", id), D1(M_D1, W_D1)
    {};

    // The reduced amplitude corresponds to the S-wave contact-like interaction
    // however it receives contributions from the propagator of the Z meson
    inline complex<double> reduced_amplitude(int i, int j)
    {
        double p_pion = _kinematics->decay_momentum_c(_s, _sab);

        // sab is assumed to be DsD channel
        complex<double> A_S = XI * (_hS / F_PION) * (1./sqrt(6.)) * sqrt(p_pion*p_pion + _mc2);
        complex<double> A_D = XI * (_hD / F_PION) * (1./sqrt(6.)) * pow( p_pion, 4.);

        // Multiply by the propagator of the D1
        A_S *= D1.eval(_sac);
        A_D *= D1.eval(_sac);

        // Can use the debug flag to turn on and off the S-wave coupling
        if (_debug == 1) A_S = 0.;

        // Being D-wave we get the appropriate projector
        // Assume pion (particle c) defined the +z direction
        return A_S * (i==j) + A_D * (3.* (i == 3)*(j ==3) - (i == j));
    };

    // -----------------------------------------------------------------------
    private:

    // We have DsDpi coupling for the S-wave and the D-wave
    double _hS = HP_S, _hD = HP_D;

    // This channel has a D1 resonance in the Ds pi channel
    relativistic_BW D1;
};

#endif