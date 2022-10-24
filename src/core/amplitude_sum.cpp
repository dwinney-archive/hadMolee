// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude_sum.hpp"

// ---------------------------------------------------------------------------
// The amplitude of a sum is simply the sum of constituent amplitudes 

complex<double> amplitude_sum::reduced_amplitude(cartesian_index i, cartesian_index j)
{
    complex<double> sum = 0.;
    for (amplitude* amp : _amps)
    {
        // Have to make sure to feed energy values to component amplitudes
        amp->update(_s, _sab, _sbc);

        sum += amp->reduced_amplitude(i, j);
    };

    return sum;
};

// ---------------------------------------------------------------------------
// Make sure all amplitudes being summed are compatible

bool amplitude_sum::check_compatibility(amplitude* amp)
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