// Wrapper class to call box loop function in LoopTools
// This needs to be compiled as its own seperate library because it needs 
// to be linked with gFortran
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef LT_BOX
#define LT_BOX

#include "box.hpp"
#include "clooptools.h"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Relativistic triangle calling LoopTools

    class LT_box : public box
    {
        // -----------------------------------------------------------------------

        public:

        LT_box()
        : box()
        {};

        LT_box(std::array<double,4> external_masses, std::array<double,4> internal_masses)
        : box(external_masses, internal_masses)
        {};

        // Evaluate by calling the B0 function
        std::complex<double> eval();
    };
};

#endif