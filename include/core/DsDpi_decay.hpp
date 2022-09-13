// Implementation of the cross-section for e+e- -> DsDpi reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "decay_process.hpp"
#include "propagator.hpp"

#ifndef DSDPI
#define DSDPI

class DsDpi_decay : public decay_process
{
    public: 
    
    // Constructor does not need anything
    // Particles assumed to be in the order 
    // a - D*
    // b - D
    // c - pi
    DsDpi_decay() 
    : decay_process(_masses, _labels, 0, "DsDpi"),
      D1_prop(M_D1, GAM_D1)
    {};


    // Output the sum of diagrams
    inline double amplitude_squared(double s, double sab, double sbc)
    {
        update(s, sab, sbc);
        return dwave_tree_squared(); 
    };

    protected:

    // ------------------------------------------------------------------------
    // D-wave contributions
    
    nonrelativistic_BW D1_prop;

    double dwave_tree_squared();

    // ------------------------------------------------------------------------
    // S-wave contributions

    // Output the modulus of pion momentum in the Y rest frame
    inline double pion_momentum()       { return sqrt(Kallen(_s, _sab, _mc2) / (4. * _s)); };
    inline double pion_momentum_squared(){ return      Kallen(_s, _sab, _mc2) / (4. * _s); };

    static const array<string,3> _labels;
    static const array<double,3> _masses;

    // Constants specific to this amplitude
    static const double GAM_D1; // (charged) D1 width in GeV
};

#endif