// Wrapper class to call triangle loop function in LoopTools
// This needs to be compiled as its own seperate library because it needs 
// to be linked with gFortran
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef LT_TRIANGLE
#define LT_TRIANGLE

#include "triangle.hpp"
#include "clooptools.h"

// ---------------------------------------------------------------------------
// Relativistic triangle calling LoopTools

class LT_triangle : public triangle
{
    // -----------------------------------------------------------------------

    public:

    LT_triangle()
    : triangle()
    {};

    LT_triangle(array<double,3> external_masses, array<double,3> internal_masses)
    : triangle(external_masses, internal_masses)
    {};

    // Evaluate by calling the C0C function
    complex<double> eval();
};

#endif