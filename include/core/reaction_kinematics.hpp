// Class to contain all relevant kinematic quantities for the process:
// e+ e- -> Y -> a b c 
// with particles a being a vector
// and b and c pseudo scalars
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef KINEMATICS
#define KINEMATICS

#include "constants.hpp"

class reaction_kinematics
{
    // -----------------------------------------------------------------------
    public:
    // Empty constructor as default
    reaction_kinematics()
    {
        set_particle_masses(0., 0., 0.);
    };

    // Specify the final state masses
    reaction_kinematics(double ma, double mb, double mc)
    {
        set_particle_masses(ma, mb, mc);
    };

    inline void set_particle_masses(double ma, double mb, double mc)
    { 
         _ma  = ma;     _mb = mb;    _mc  = mc;
         _ma2 = ma*ma; _mb2 = mb*mb; _mc2 = mc*mc;
    };

    // Set particle labels (useful for plotting)
    inline void set_particle_labels(string a, string b, string c){ _a = a; _b = b; _c = c; };

    // -----------------------------------------------------------------------
    private:
    
    // Particles can carry a string label
    string _a = "a", _b = "b", _c = "c";

    // Overall center-of-mass energy
    double _W, _s;

    // Sub-channel invariant masses
    double _sab, _sbc, _sac;

    // Masses of the particles
    double _ma,  _mb,  _mc;
    double _ma2, _mb2, _mc2;
};

#endif