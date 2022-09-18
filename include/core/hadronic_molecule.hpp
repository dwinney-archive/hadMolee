// Container class holding all quantities relevant for a hadronic molecule state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef MOLECULE
#define MOLECULE

#include "constants.hpp"
#include "propagator.hpp"

class hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    hadronic_molecule(string id = "exotic")
    : _id(id)
    {};

    virtual complex<double> propagator(double s){ return 1.; };

    // -----------------------------------------------------------------------
    private:

    // Name identifier
    string _id;

    // Masses
    double _bare_mass, _renormalized_mass;

    // Constituents masses 
    double _m1, _m2;
    double _molecule_coupling;
};

#endif