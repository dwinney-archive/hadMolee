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

#include "Math/Functor.h"
#include "Math/Derivator.h"
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

    // Evaluate the propagator
    // For standardization the input argument is always asssumed to be s, take the square root internally if we need E
    virtual complex<double> propagator(double s){ return 1.; };

    // String identifier
    inline string get_id(){ return _id; };

    // Output the saved coupling to the constituent channel
    virtual inline double coupling(){ return _bare_coupling; };

    // Setting utilitites for diffferent parameters that can float
    virtual inline void set_pole_mass(double x){ _renormalized_mass = x; };
    virtual inline void set_coupling( double x){ _bare_coupling     = x; };

    // Set pole mass ,coupling, and non-mol width in a single call
    virtual inline void set_parameters(array<double,3> pars)
    {
        set_pole_mass(pars[0]);
        set_coupling( pars[1]);
        _nonmolecular_width = pars[3];
    };  

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
    double _bare_coupling = 1.;
    inline double reduced_mass(){ return _m1 * _m2 / (_m1 + _m2); };
    inline double mass_difference(double E){ return E - _m1 - _m2; };

    // Intermediate state threshold openings of the constituent channel
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
        _bare_mass          = M_ZC3900;
        _total_width        = W_ZC3900;
        
        // Coupling taken from [1]
        _bare_coupling      = C_Z;
        
        // residual width taken to recover the full PDG width at the pole
        _nonmolecular_width = _total_width - 2.* imag(self_energy(_bare_mass));
    };

    // The propagator gains contributions from the self-energy
    complex<double> propagator(double s);

    // Self-energy from bubble diagram of D* D scattering 
    complex<double> self_energy(double E);

    // -----------------------------------------------------------------------
    private:

    // Total width of the Z from the PDG
    double _total_width;
};

// ---------------------------------------------------------------------------
// Implementation of D1 D molecule for the Y(4260)

class D1D_molecule : public hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    D1D_molecule(string id = "Y(4260)")
    : hadronic_molecule(M_D1, M_D, id)
    {
        // Mass and Width from PDG
        _bare_mass          = M_Y4260;
        
        // Coupling taken from [1]
        _bare_coupling      = C_Y;

        // Set up the derivator 
        wsigma = ROOT::Math::Functor1D(this, &D1D_molecule::resigma);
        dsigma.SetFunction(wsigma);
    };

    // The propagator gains contributions from the self-energy
    complex<double> propagator(double s);

    // When we change the pole mass, we must recalculate renormalization quantities
    void set_pole_mass(double x)
    {
        _renormalized_mass = x;

        // Precalculate relevant quantities
        _reS  = resigma(_renormalized_mass);
        _redS = dsigma.Eval(_renormalized_mass);
        _Z    = 1. / (1. - _redS);
    };

    private:

    // Self-energy from bubble of D1 D scattering and dressed with elastic scattering
    // renomalized
    complex<double> self_energy(double E);

    // Bare self-energy just from the bubble of D1 D scattering
    complex<double> sigma(double E);
    inline double resigma(double E){ return real(sigma(E)); };

    // Need to be able to calculate the derivative of the above self-energy
    ROOT::Math::Functor1D wsigma;
    ROOT::Math::Derivator dsigma;

    double _Z, _reS, _redS;
};

#endif