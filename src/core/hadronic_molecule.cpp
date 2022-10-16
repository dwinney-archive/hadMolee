// Container class holding all quantities relevant for a hadronic molecule state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "hadronic_molecule.hpp"

// ---------------------------------------------------------------------------
// Z MESON

// Z meson self-energy from D* D scattering,
// this is the "bare" self-energy with no rescattering contribution for D*D elastic diagram
complex<double> DsD_molecule::self_energy(double E)
{
    double eps = mass_difference(E);
    double mu  = reduced_mass();
    double z   = _bare_coupling;
    
    return (z*z / (8.*PI)) * sqrt(2.*mu*mu*mu*abs(eps)) * ( XI*(eps>=0) - XR*(eps<0) );
};

// Propagator recieves contributions from the self-energy above and the constant width
complex<double> DsD_molecule::propagator(double E)
{
    complex<double> D = E - _bare_mass + self_energy(E) + XI*_nonmolecular_width/2.;
    
    return XI / (2.*D);
};  

// ---------------------------------------------------------------------------
// Y MESON

// Bare self-energy, i.e. the simple D1 D bubble in MSbar renormalization
// Because the D1 has a width, we include this explicitly as a constant
// This renders the self-energy smoothly defined for all E
complex<double> D1D_molecule::sigma(double E)
{
    double eps = mass_difference(E);
    double mu  = reduced_mass();
    double y   = _bare_coupling;

    return (-XI*y*y / (8.*PI)) * sqrt( 2.*mu*mu*mu* (eps + XI*W_D1/2.) );
};

// Calcualte the renormalized self-energy by subtracting away the real parts of the function and its first derivative
complex<double> D1D_molecule::self_energy(double E)
{
    return sigma(E) - _reS - (E - _renormalized_mass) * _redS;
};

// Full propagator
complex<double> D1D_molecule::propagator(double E)
{
    complex<double> D = E - _renormalized_mass - _Z*self_energy(E) + XI*_nonmolecular_width/2.;
    return XI * _Z / (2. * D);
};
