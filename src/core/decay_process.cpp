// Abstract class for an decay_process of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "decay_process.hpp"


// ---------------------------------------------------------------------------
// Doubly differential partial-width
double decay_process::d2Gamma(double s, double sab, double sbc)
{
    if ( !in_physical_region(s, sab, sbc) ) 
    {
        warning(get_id(), "Evaluating decay_process outside physical region! Returning 0...");
        return 0.;
    };

    // Update the saved energy values for easier access later
    update(s, sab, sbc);

    double amp_squared = amplitude_squared(s, sab, sbc); 

    // General prefactors for 1->3 decay width in GeV
    double prefactors = 1. / (32.*pow(2.*PI*sqrt(s), 3.));

    // Additional prefactors since we incorporate the leptonic production here 
    // i.e. the e+ e- -> \gamma vertex at high energies gives
    prefactors /= E*E * pow(sqrt(s), 3.);

    return amp_squared * prefactors;
};

// ---------------------------------------------------------------------------
// Singly differentrial partial-widths

// In the ab subsystem
double decay_process::dGamma_ab(double s, double sab)
{
    auto F = [&](double sbc)
    {
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = sbc_from_sab(s, sab, -1.);
    double sbc_max = sbc_from_sab(s, sab, +1.);

    return ig.Integral(sbc_min, sbc_max);
};

// In the bc subsystem
double decay_process::dGamma_bc(double s, double sbc)
{
    auto F = [&](double sab)
    {
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = sab_from_sbc(s, sbc, -1.);
    double sab_max = sab_from_sbc(s, sbc, +1.);

    return ig.Integral(sab_min, sab_max);
};

// In the ac subsystem
double decay_process::dGamma_ac(double s, double sac)
{
    auto F = [&](double sbc)
    {
        double sac = _ma2 + _mb2 + _mc2 + s - sac - sbc;
        return d2Gamma(s, sac, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = sbc_from_sac(s, sac, -1.);
    double sbc_max = sbc_from_sac(s, sac, +1.);

    return ig.Integral(sbc_min, sbc_max);
};

// ---------------------------------------------------------------------------
// Fully integrated decay width
double decay_process::Gamma(double s)
{
    auto F = [&](double sab)
    {
        return dGamma_ab(s, sab);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = (_ma + _mb)*(_ma + _mb);
    double sab_max = (  s - _mc)*(  s - _mc);

    return ig.Integral(sab_min, sab_max);
}