// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude_sum.hpp"

// ---------------------------------------------------------------------------
// Check if two amplitudes (in smart_ptr form) are compatible to be added 
bool hadMolee::are_compatible(std::shared_ptr<amplitude> a, std::shared_ptr<amplitude> b)
{
    // Compatible amplitudes have the same kineamtics and Y meson pointers
    if (a->_kinematics != b->_kinematics)
    {
        warning("amplitude_sum", 
                "Amplitudes " + a->get_id() + " and " + b->get_id() + " contain different stored kinematics instances (" + a->_kinematics->get_id() + " and " + b->_kinematics->get_id() + "). Returning null...");
        return false;
    } 
    if (a->_V != b->_V) 
    {
        warning("amplitude_sum", 
                "Amplitudes " + a->get_id() + " and " + b->get_id() + " contain different e+ e- line-shape instances (" + a->_V->get_id() + " and " + b->_V->get_id() + "). Skipping null...");
        return false;
    }
    
    return true;
};

// ---------------------------------------------------------------------------
// The amplitude of a sum is simply the sum of constituent amplitudes 

std::complex<double> hadMolee::amplitude_sum::reduced_amplitude(cartesian_index i, cartesian_index j)
{
    std::complex<double> sum = 0.;
    for (auto amp : _amps)
    {
        // Have to make sure to feed energy values to component amplitudes
        amp->update(_s, _sab, _sbc);

        sum += amp->reduced_amplitude(i, j);
    };

    return sum;
};

// ---------------------------------------------------------------------------
// Make sure all amplitudes being summed are compatible
// This is efffectively the same as the "are_compatible" except comparisions are made against the interally stored instances 

bool hadMolee::amplitude_sum::is_compatible(std::shared_ptr<amplitude> amp)
{
    // Compatible amplitudes have the same kineamtics and Y meson pointers
    if (amp->_kinematics != _kinematics)
    {
        warning("amplitude_sum", 
                "New amplitude " + amp->get_id() + " does not contain the stored kinematics instance (" + _kinematics->get_id() + "). \n Skipping amplitude...");
        return false;
    } 
    if (amp->_V != _V) 
    {
        warning("amplitude_sum", 
                "New amplitude " + amp->get_id() + " does not contain the stored hadronic_molecule instance (" + _V->get_id() + "). \n Skipping amplitude...");
        return false;
    }
    
    return true;
};