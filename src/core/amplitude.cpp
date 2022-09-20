// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// At each energy step we cach all the components of the reduced amplitude tensor to minimize re-calculation

inline void amplitude::check_cache()
{
    bool need_recalculate;
    need_recalculate = (_cached_amplitudes.size() != 3 || _cached_amplitudes[0].size() != 3 ) 
                        || (abs(_cached_s - _s) > _cache_tolerance) 
                        || (abs(_cached_sab - _sab) > _cache_tolerance) 
                        || (abs(_cached_sbc - _sbc) > _cache_tolerance);
    if (need_recalculate)
    {
        _cached_amplitudes.clear();
        for (int i = 1; i < 4; i++)
        {
            vector<complex<double>> i_slice;
            for (int j = 1; j < 4; j++)
            {
                i_slice.push_back( reduced_amplitude(i, j) );
            }

            _cached_amplitudes.push_back(i_slice);
        };
    };
};

// ---------------------------------------------------------------------------
// Spin-summed amplitude squared. 
double amplitude::probability_distribution(double s, double sab, double sbc)
{
    // Update the saved energy values for easier access later
    update(s, sab, sbc);

    // Check if we need to update our precalculated amplitudes 
    check_cache();

    double sum = 0;

    // Contract over Cartesian indices
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                complex<double> x;
                x  = _kinematics->production_tensor(i+1, j+1, s);
                x *= _Y->propagator(s);
                x *=      _cached_amplitudes[j][k];
                x *= conj(_cached_amplitudes[k][i]);

                if (imag(x) > 1.E-5) warning("probability_distribution", "Spin sum squared is imaginary!!!!");

                sum += real(x);
            };
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
        warning(get_id(), "Evaluating amplitude outside physical region! Returning 0...");
        return 0.;
    };

    double amp_squared = probability_distribution(s, sab, sbc); 

    // Average over the spins of the inital electron-positron pair
    amp_squared /= 4.;

    // General prefactors for 1->3 decay width in GeV
    double prefactors = 1. / (32.*pow(2.*PI*sqrt(s), 3.));

    if (_normalize) prefactors *= _normalization;

    return amp_squared * prefactors;
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

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = _kinematics->sbc_from_sab(s, sab, -1.);
    double sbc_max = _kinematics->sbc_from_sab(s, sab, +1.);

    if (abs(sbc_min - sbc_max) < EPS) return 0.;
    return ig.Integral(sbc_min, sbc_max);
};

// In the bc subsystem
double amplitude::dGamma_bc(double s, double sbc)
{
    auto F = [&](double sab)
    {
        return d2Gamma(s, sab, sbc);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = _kinematics->sab_from_sbc(s, sbc, -1.);
    double sab_max = _kinematics->sab_from_sbc(s, sbc, +1.);

    if (abs(sab_min - sab_max) < EPS) return 0.;
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

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sbc_min = _kinematics->sbc_from_sac(s, sac, -1.);
    double sbc_max = _kinematics->sbc_from_sac(s, sac, +1.);

    if (abs(sbc_min - sbc_max) < EPS) return 0.;
    return ig.Integral(sbc_min, sbc_max);
};

// ---------------------------------------------------------------------------
// Fully integrated decay width
double amplitude::Gamma(double s)
{
    auto F = [&](double sab)
    {
        return dGamma_ab(s, sab);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS31);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double sab_min = (_ma + _mb)*(_ma + _mb);
    double sab_max = ( sqrt(s) - _mc)*( sqrt(s) - _mc);

    return ig.Integral(sab_min, sab_max);
}