// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// At each energy step we cach all the components of the reduced amplitude tensor to minimize re-calculation

inline void amplitude::check_decay_cache()
{
    bool need_recalculate;
    need_recalculate =     (abs(_cached_s - _s) > _cache_tolerance) 
                        || (abs(_cached_sab - _sab) > _cache_tolerance) 
                        || (abs(_cached_sbc - _sbc) > _cache_tolerance);
    if (need_recalculate)
    {
        // We neeed to calcualte the amplitude squared summed over the final state polarizations
        // First lets calcualte and save the reduced amplitude for i,j=0,1,2
        array<array<complex<double>,3>,3> amp;

        for (auto i : C_INDICES)
        {
            for (auto j : C_INDICES)
            {
                amp[i][j] = reduced_amplitude(i, j);
            }
        };

        //Now we can square and sum without having to calculate twice
        for (auto i : C_INDICES)
        {
            for (auto j : C_INDICES)
            {
                // K is the external polarization
                for (auto k : C_INDICES)
                {
                    complex<double> x;
                    x  =      amp[i][k];
                    // Polarization sum over the final state vector gives a delta function fixing k here
                    x *= conj(amp[j][k]);
                    if (!is_zero( imag(x) )) warning("check_decay_cache", "Reduced amplitude squared is imaginary!");
                    _cached_decay_tensor[i][j] = real(x);
                };
            }
        };
    };
};

// ---------------------------------------------------------------------------
// Spin-summed amplitude squared. 

// This is split into two pieces depending on whether or not to include the e+e- components
// The decay distribution is the amplitude squared for the process V->abc
double amplitude::decay_distribution(double s, double sab, double sbc)
{
    // Update the saved energy values for easier access later
    update(s, sab, sbc);

    // Contract over Cartesian indices
    double sum = 0;
    for (auto i : C_INDICES)
    {
        for (auto j : C_INDICES)
        {
            int V_polarization = delta(i, j) /* - delta(i,z)*delta(j,z) */;
            if (V_polarization == 0) continue;

            // V -> abc 
            double Asqr = _cached_decay_tensor[i][j];
            if ( is_zero(Asqr) ) continue;

            sum += V_polarization * Asqr;
        };
    };

    return sum;
};


// The scattering distribution is the amplitude squared for the process e+e- ->abc
double amplitude::scattering_distribution(double s, double sab, double sbc)
{
    // Update the saved energy values for easier access later
    update(s, sab, sbc);

    // Contract over Cartesian indices
    double sum = 0;
    for (auto i : C_INDICES)
    {
        for (auto j : C_INDICES)
        {
            int production = delta(i, j) - delta(i,z)*delta(j,z);
            if (production == 0) continue;

            complex<double> x = 1.;
            //e+ e- -> gamma 
            // x *= _kinematics->ee_to_gamma(s) * production;
            x = XR * production;

            // gamma -> V
            // x *= norm(_V->photon_coupling() * _V->propagator(s));

            // V -> abc 
            x *= _cached_decay_tensor[i][j];

            if (!is_zero( imag(x) )) warning("probability_distribution", "Amplitude squared is imaginary!");
            sum += real(x);
        };
    };

    return sum;
};


// ---------------------------------------------------------------------------
// Doubly differential partial-width
double amplitude::d2Gamma(double s, double sab, double sbc)
{
    if ( !_kinematics->in_physical_region(s, sab, sbc) ) 
    {
        // warning(get_id(), "Evaluating amplitude outside physical region! Returning 0...");
        return 0.;
    };

    double amp_squared = decay_distribution(s, sab, sbc); 

    // // Average over the spins of the inital electron-positron pair
    amp_squared /= 2.;

    // General prefactors for 1->3 decay width in GeV
    double prefactors = 1. / (32.*pow(2.*PI*sqrt(s), 3.));

    if (_normalize) prefactors *= _normalization;

    return amp_squared * prefactors * 1.E3; // In MeV
};

// ---------------------------------------------------------------------------
// Singly differentrial partial-widths

// In the ab subsystem
double amplitude::dGamma_ab(double s, double sab)
{
    auto F = [&](double sbc)
    {
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = _kinematics->sbc_from_sab(s, sab, -1.);
    double sbc_max = _kinematics->sbc_from_sab(s, sab, +1.);

    if ( is_equal(sbc_min, sbc_max) ) return 0.;
    return ig.Integral(sbc_min, sbc_max);
};

// In the bc subsystem
double amplitude::dGamma_bc(double s, double sbc)
{
    auto F = [&](double sab)
    {
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = _kinematics->sab_from_sbc(s, sbc, -1.);
    double sab_max = _kinematics->sab_from_sbc(s, sbc, +1.);

    if ( is_equal(sab_min, sab_max) ) return 0.;
    return ig.Integral(sab_min, sab_max);
};

// In the ac subsystem
double amplitude::dGamma_ac(double s, double sac)
{
    auto F = [&](double sbc)
    {
        double sab = _ma2 + _mb2 + _mc2 + s - sac - sbc;
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = _kinematics->sbc_from_sac(s, sac, -1.);
    double sbc_max = _kinematics->sbc_from_sac(s, sac, +1.);

    if ( is_equal(sbc_min, sbc_max) ) return 0.;
    return ig.Integral(sbc_min, sbc_max);
};

// Alias for three different subchannels, specify with an argument
// Third argument is expected to be the correct subchannel energy
double amplitude::dGamma(subchannel chan, double s, double sigma)
{
    if ( !_kinematics->in_physical_region(s, sigma, chan) ) 
    {
        return 0.;
    };

    switch (chan)
    {
        case ab: return dGamma_ab(s, sigma);
        case bc: return dGamma_bc(s, sigma);
        case ac: return dGamma_ac(s, sigma);
        default: return 0.;
    };

    return 0.;
};

// ---------------------------------------------------------------------------
// Fully integrated decay width
double amplitude::Gamma(double s)
{
    auto F = [&](double sab)
    {
        return dGamma_ab(s, sab);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = (_ma + _mb)*(_ma + _mb);
    double sab_max = ( sqrt(s) - _mc)*( sqrt(s) - _mc);

    return ig.Integral(sab_min, sab_max);
}