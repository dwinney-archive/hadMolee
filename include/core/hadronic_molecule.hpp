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
#include "charmoniumlike.hpp"

// ---------------------------------------------------------------------------
// Generic class, specific implementations for the Y and Z mesons are given below

class hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    // Constructor requires setting the masses of constituent
    // Set the constituent channel masses
    hadronic_molecule(double m1, double m2, int npars)
    : _m1(m1), _m2(m2), _npars(npars)
    {};

    // Evaluate the propagator
    // For standardization the input argument is always asssumed to be s, take the square root internally if we need E
    virtual complex<double> propagator(double s){ return 1.; };

    // Output the saved coupling to the constituent channel
    virtual inline double molecular_coupling(){ return _coupling; };

    // // Setting utilitites for diffferent parameters that can float
    // virtual inline void set_pole_mass(double x){ _renormalized_mass = x; };
    // virtual inline void set_coupling( double x){ _bare_coupling     = x; };

    // Set pole mass ,coupling, and non-mol width in a single call
    virtual inline void set_parameters(vector<double> pars)
    {
        check_size(pars);
        return;
    };  

    // -----------------------------------------------------------------------
    protected:

    // Name identifier
    string _id;

    // Masses of the molecule
    double _pole_mass;

    // Width coming from decays other than m1 m2 final state
    double _nonmol_width;

    // Constituents masses 
    double _m1, _m2;
    double _coupling = 1.;
    inline double reduced_mass(){ return _m1 * _m2 / (_m1 + _m2); };
    inline double mass_difference(double E){ return E - _m1 - _m2; };

    // Intermediate state threshold openings of the constituent channel
    inline double Wth(){ return _m1 + _m2; };
    inline double sth(){ return Wth()*Wth(); };

    // Self-energy loop function
    virtual complex<double> self_energy(double x){ return 0.; };

    void check_size(vector<double> pars)
    {
        if (pars.size() != _npars)
        {
            warning("hadronic_molecule", "Wrong number of parameters given! Expected 3 but recieved " + to_string(pars.size()) + ". Results may vary...");
        };
    }
    int _npars = 0;
};

// ---------------------------------------------------------------------------
// D* D molecule relevant for the Z meson

class DsD_molecule : public hadronic_molecule
{
    // -----------------------------------------------------------------------
    public:

    DsD_molecule()
    : hadronic_molecule(M_DSTAR, M_D, 0)
    {
        // Mass and Width from PDG
        _pole_mass          = M_ZC3900;
        _total_width        = W_ZC3900;
        
        // Coupling taken from [1]
        _coupling      = ZBARE_QQ2016;
        
        // residual width taken to recover the full PDG width at the pole
        _nonmol_width = _total_width - 2.* imag(self_energy(_pole_mass));
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

class D1D_molecule : public hadronic_molecule, public charmoniumlike
{
    // -----------------------------------------------------------------------
    public:

    D1D_molecule(string id = "Y(4260)")
    : hadronic_molecule(M_D1, M_D, 4), charmoniumlike(id)
    {
        // Set up the derivator 
        wsigma = ROOT::Math::Functor1D(this, &D1D_molecule::resigma);
        dsigma.SetFunction(wsigma);

        // Default parameters use the QQ2016 values
        set_parameters({MY_QQ2016, YBARE_QQ2016, YNM_WIDTH_QQ2016, F_Y_QQ2016});
    };

    // The propagator gains contributions from the self-energy
    complex<double> propagator(double s);

    // Since Y-meson is also charmonium-like it requires a photon coupling
    complex<double> photon_coupling()
    {
        return XI * E * _pole_mass*_pole_mass / _fY;
    };

    // When we change the pole mass, we must recalculate renormalization quantities
    // Set pole mass ,coupling, and non-mol width in a single call
    inline void set_parameters(vector<double> pars)
    {
        check_size(pars);

        _pole_mass    = pars[0];
        _coupling     = pars[1];
        _nonmol_width = pars[2];
        _fY           = pars[3];

        // Precalculate relevant quantities
        _reS  = resigma(_pole_mass);
        _redS = dsigma.Eval(_pole_mass);
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

    double _Z, _reS, _redS, _fY;
};

#endif