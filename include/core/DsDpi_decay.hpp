// Implementation of the cross-section for e+e- -> DsDpi reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "decay_process.hpp"

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
    : decay_process(_masses, _labels, 0, "DsDpi")
    {};


    // Output the sum of diagrams
    inline double amplitude_squared(double s, double sab, double sbc){ return 1.; };
    protected:

    static const array<string,3> _labels;
    static const array<double,3> _masses;
};

#endif