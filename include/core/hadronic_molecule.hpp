// Container class holding all quantities relevant for a hadronic molecule state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// References:

// [1] arXiv:1310.2190 [hep-ph]
// ---------------------------------------------------------------------------


#ifndef MOLECULE
#define MOLECULE

#include "constants.hpp"

// ---------------------------------------------------------------------------
// Generic class, specific implementations for the Y and Z mesons are given below

class hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    // Empty constructor with nothign set except maybe an id
    hadronic_molecule(string id = "exotic")
    {};

    // Set the constituent channel masses
    hadronic_molecule(double m1, double m2, string id = "exotic")
    : _m1(m1), _m2(m2), _id(id)
    {};

    // Evaluate the propagator, the variable here could be E or s depending on if relativsitic or not
    virtual complex<double> propagator(double x){ return 1.; };

    // -----------------------------------------------------------------------
    protected:

    // Name identifier
    string _id;

    // Masses of the molecule
    double _bare_mass, _renormalized_mass;
    double _renomalization_constant;

    // Width coming from decays other than m1 m2 final state
    double _nonmolecular_width;

    // Constituents masses 
    double _m1, _m2;
    double _bare_coupling;
    inline double reduced_mass(){ return _m1 * _m2 / (_m1 + _m2); };
    inline double mass_difference(double E){ return E - _m1 - _m2; };

    // Threshold openings of the constituent channel
    inline double Wth(){ return _m1 + _m2; };
    inline double sth(){ return Wth()*Wth(); };

    // Self-energy loop function
    virtual complex<double> self_energy(double x){ return 0.; };
};

// ---------------------------------------------------------------------------
// D* D molecule relevant for the Z meson

class DsD_molecule : public hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    DsD_molecule(string id = "Zc(3900)")
    : hadronic_molecule(M_DSTAR, M_D, id)
    {
        // Mass and Width from PDG
        _bare_mass          = M_Z;
        _total_width        = 28.4E-3;
        
        // Coupling taken from [1]
        _bare_coupling      = 0.77;
        
        // residual width taken to recover the full PDG width at the pole
        _nonmolecular_width = _total_width - 2.* imag(self_energy(_bare_mass));
    };

    // The propagator gains contributions from the self-energy
    complex<double> propagator(double E);

    // Self-energy from bubble diagram of D* D scattering 
    complex<double> self_energy(double E);

    // -----------------------------------------------------------------------
    private:

    // Total width of the Z from the PDG
    double _total_width;

};

#endif