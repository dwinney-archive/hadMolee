// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Make sure all amplitudes being summed are compatible
    // This is efffectively the same as the "are_compatible" except comparisions are made against the interally stored instances 

    bool are_compatible(amplitude a, amplitude b)
    {
        // Compatible amplitudes have the same kineamtics and Y meson pointers
        if (a->_kinematics != b->_kinematics)
        {
            warning("amplitude", 
                    "Amplitudes " + a->get_id() + " and " + b->get_id() + " contain different stored kinematics instances (" + a->_kinematics->get_id() + " and " + b->_kinematics->get_id() + "). Returning null...");
            return false;
        } 
        if (a->_V != b->_V) 
        {
            warning("amplitude", 
                    "Amplitudes " + a->get_id() + " and " + b->get_id() + " contain different e+ e- line-shape instances (" + a->_V->get_id() + " and " + b->_V->get_id() + "). Skipping null...");
            return false;
        }
        
        return true;
    };

    // Smart pointer "constructor" for a sum of amplitudes 
    amplitude operator+(amplitude a, amplitude b)
    {
        if ( !are_compatible(a, b) ) return nullptr;

        auto sum = std::make_shared<amplitude_base>(amplitude_key(), a->_kinematics, a->_V, a->get_id() + " + " + b->get_id());

        // If the constituent amplitudes are already sums, add the vector of contituents
        // if theyre a 'bare' amplitude add the amplitude itself
        (a->is_sum()) ? ( sum->add({a->_sub_amps}) ) : ( sum->add(a) );
        (b->is_sum()) ? ( sum->add({b->_sub_amps}) ) : ( sum->add(b) );

        return sum;
    };

    void operator+=(amplitude a, amplitude b)
    {
        if (!a->is_sum())
        {
            warning("amplitude", "Tried adding to an amplitude (" + a->get_id() + ") which is not already a sum!" + 
                                "\nInitialize a sum first by using the binary + operator then increment with += to avoid unexpected behavior." + 
                                "\nReturning without change...");
            return;
        };

        a->add(b);
    };

    // ---------------------------------------------------------------------------
    // These are for amplitude summing from existing amplitudes

    bool amplitude_base::is_compatible(amplitude amp)
    {
        // Compatible amplitudes have the same kineamtics and Y meson pointers
        if (amp->_kinematics != _kinematics)
        {
            warning("amplitude", 
                    "New amplitude " + amp->get_id() + " does not contain the stored kinematics instance (" + _kinematics->get_id() + "). \n Skipping amplitude...");
            return false;
        } 
        if (amp->_V != _V) 
        {
            warning("amplitude", 
                    "New amplitude " + amp->get_id() + " does not contain the stored hadronic_molecule instance (" + _V->get_id() + "). \n Skipping amplitude...");
            return false;
        }
        
        return true;
    };

    // ---------------------------------------------------------------------------
    // At each energy step we cach all the components of the reduced amplitude tensor to minimize re-calculation

    void amplitude_base::check_decay_cache()
    {
        bool need_recalculate;
        need_recalculate =     (std::abs(_cached_s   - _s   ) > _cache_tolerance) 
                            || (std::abs(_cached_sab - _sab ) > _cache_tolerance) 
                            || (std::abs(_cached_sbc - _sbc ) > _cache_tolerance)
                            || (std::abs(_cached_cos - _cos ) > _cache_tolerance);
        if (need_recalculate)
        {
            // We neeed to calcualte the amplitude squared summed over the final state polarizations
            // First lets calcualte and save the reduced amplitude for i,j=0,1,2
            std::array<std::array<complex,3>,3> amp;

            for (auto i : C_INDICES)
            {
                for (auto j : C_INDICES)
                {
                    amp[i][j] = reduced_amplitude_checked(i, j);
                }
            };

            //Now we can square and sum without having to calculate twice
            for (auto i : C_INDICES)
            {
                for (auto j : C_INDICES)
                {
                    complex x = 0.;

                    // K is the external polarization
                    for (auto k : C_INDICES)
                    {
                        // Polarization sum over the final state vector gives a delta function fixing k here
                        x += amp[i][k] * conj(amp[j][k]);
                    };

                    if (!is_zero( imag(x) ))
                    {
                        warning("check_decay_cache", "Reduced amplitude squared is imaginary!");
                    };
                    _cached_decay_tensor[i][j] = real(x);
                }
            };
        };
    };

    // Allocate an aggragated vector of parameters to individual amplitudes
    void amplitude_base::set_parameters(std::vector<double> x)
    {
        check_nParams(x);

        // If called from empty amplitude, do nothing
        if (!is_sum()) return;

        // Else allocate parameters to the individual amplitudes
        int N = 0;
        for (auto amp : _sub_amps)
        {
            // Extract the subvector corresponding to the i-th amplitude
            auto start = x.begin() + N;
            auto end   = x.begin() + N + amp->N_parameters();
            std::vector<double> pars(start, end);

            amp->set_parameters(pars);
            N += amp->N_parameters();
        };

        // At the end check that the number of params allocated is the same as expected
        // However is this happpens the damage is done so we just send out a warning...
        if (N != N_parameters())
        {
            std::cout << "Warning! amplitude::set_params() : Number of parameters allocated doesnt match those expected..." << std::endl;
        };
    };

    // ---------------------------------------------------------------------------
    // The amplitude of a sum is simply the sum of constituent amplitudes 

    complex amplitude_base::reduced_amplitude(cartesian_index i, cartesian_index j)
    {
        // IF called without any store amplitudes, throw error
        if (!is_sum()) return std::nan("");
        
        // Else sum all constituent amplitudes together
        complex sum = 0.;
        for (auto amp : _sub_amps)
        {
            // Have to make sure to feed energy values to component amplitudes
            amp->update(_s, _sab, _sbc, _cos);
            sum += amp->reduced_amplitude_checked(i, j);
        };
        return sum;
    };

    // ---------------------------------------------------------------------------
    // Spin-summed amplitude squared. 

    // The scattering distribution is the amplitude squared for the process e+e- ->abc
    double amplitude_base::decay_distribution(double s, double sab, double sbc, double cos)
    {
        // Update the saved energy values for easier access later
        update(s, sab, sbc, cos);

        // Contract over Cartesian indices
        double sum = 0;
        for (auto i : C_INDICES)
        {
            for (auto j : C_INDICES)
            {
                complex x; //Output

                // --------------------------
                // e+ e- -> gamma component

                // Production tensor of the e+ e- polarization component
                int production_pol = delta(i, j) - delta(i,z)*delta(j,z);
                if (production_pol == 0) continue;

                // ( e^2 * 2k^2 * polarization)
                double k = sqrt(s)/2; // e+ momentum in CM frame
                x  = E*E * 2*k*k * production_pol;

                // --------------------------
                // gamma -> V component

                // ( proton propagator * Vgamma coupling * V propagator)
                complex Gprime = (1/s) * _V->photon_coupling() * _V->propagator(s);
                x *= norm( Gprime );

                // --------------------------
                // Decay component

                // V -> abc 
                x *= _cached_decay_tensor[i][j];

                if (!is_zero( imag(x) )) warning("probability_distribution", "Amplitude squared is imaginary!");
                sum += real(x);
            };
        };

        return sum;
    };

    // Without specifying an angle assumes we integrate over cos
    double amplitude_base::decay_distribution(double s, double sab, double sbc)
    {
        auto dF = [&](double cos)
        {
            return decay_distribution(s, sab, sbc, cos);
        };

        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF(dF);
        ig.SetFunction(wF);

        return ig.Integral(-1., 1.);
    };

    // ---------------------------------------------------------------------------
    // Doubly differential partial-width

    double amplitude_base::d2Gamma(double s, double sab, double sbc)
    {
        if ( !_kinematics->in_physical_region(s, sab, sbc) ) 
        {
            return 0.;
        };

        // General prefactors for 1->3 decay width in GeV
        double prefactors = 1. / (32.*pow(2.*PI*sqrt(s), 3.));

        if (_normalize) prefactors *= _normalization;

        return decay_distribution(s, sab, sbc) * prefactors * 1.E3; // In MeV
    };

    // ---------------------------------------------------------------------------
    // Singly differentrial partial-widths

    // In the ab subsystem
    double amplitude_base::dGamma_ab(double s, double sab)
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
    double amplitude_base::dGamma_bc(double s, double sbc)
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
    double amplitude_base::dGamma_ac(double s, double sac)
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
    double amplitude_base::dGamma(subchannel chan, double s, double sigma)
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
    double amplitude_base::Gamma(double s)
    {
        auto F = [&](double sab)
        {
            return dGamma(ab, s, sab);
        };

        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF(F);
        ig.SetFunction(wF);

        double sab_min = (_ma + _mb)*(_ma + _mb);
        double sab_max = ( sqrt(s) - _mc)*( sqrt(s) - _mc);

        return ig.Integral(sab_min, sab_max);
    }
};