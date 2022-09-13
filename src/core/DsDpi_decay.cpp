// Implementation of the cross-section for e+e- -> DsDpi reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "DsDpi_decay.hpp"

const array<string,3> DsDpi_decay::_labels = {"#it{D}*", "#it{D}", "#pi" };
const array<double,3> DsDpi_decay::_masses = {  M_DSTAR, M_D, M_PION };

const double DsDpi_decay::GAM_D1 = 31.3E-3;

// The squared tree-level transition via D1 resonance in D-wave
double DsDpi_decay::dwave_tree_squared()
{
    // D-wave prefactor times BW for D1
    double p_pi4 = pion_momentum_squared() * pion_momentum_squared();
    return 4. * p_pi4 * D1_prop.squared(sqrt(_sac));
};