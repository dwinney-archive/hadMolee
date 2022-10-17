// Container class holding all quantities relevant for a charmoniumlike particle
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// References:

// [1] arXiv:1310.2190 [hep-ph]
// ---------------------------------------------------------------------------


#ifndef CHARMONIUMLIKE
#define CHARMONIUMLIKE

#include "constants.hpp"

class charmoniumlike
{
    // -----------------------------------------------------------------------
    public:

    // Empty constructor
    charmoniumlike(string id = "charmonium")
    {};

    virtual complex<double> propagator(double s){ return 1.; };
    virtual complex<double> photon_coupling(){    return 1.; };

    // String identifier
    inline string get_id(){ return _id; };

    // -----------------------------------------------------------------------
    protected:
    
    // Name identifier
    string _id;
};

#endif