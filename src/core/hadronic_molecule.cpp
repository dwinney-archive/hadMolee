// Container class holding all quantities relevant for a hadronic molecule state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "hadronic_molecule.hpp"

// ---------------------------------------------------------------------------
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
