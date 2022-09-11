// Class to contain all relevant kinematic quantities for the process:
// e+ e- -> Y -> a b c 
// with particles a being a vector
// and b and c pseudo scalars
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef KINEMATICS
#define KINEMATICS

#include "constants.hpp"
#include "helicities.hpp"

// Kallen triangle function
template <typename T>
inline T Kallen(T x, T y, T z)
{
    return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
};

class reaction_kinematics
{
    // -----------------------------------------------------------------------
    public:
    // Empty constructor as default
    reaction_kinematics()
    {
        set_particle_masses(0., 0., 0.);
    };

    // Specify the final state masses
    reaction_kinematics(double ma, double mb, double mc)
    {
        set_particle_masses(ma, mb, mc);
    };

    inline void set_particle_masses(double ma, double mb, double mc)
    { 
         _ma  = ma;     _mb = mb;    _mc  = mc;
         _ma2 = ma*ma; _mb2 = mb*mb; _mc2 = mc*mc;
    };
    inline array<double,3> get_masses(){ return {_ma, _mb, _mc}; };


    // Set particle labels (useful for plotting)
    inline void set_particle_labels(string a, string b, string c){ _a = a; _b = b; _c = c; };
    inline string particle_a(){ return _a; };
    inline string particle_b(){ return _b; };
    inline string particle_c(){ return _c; };


    // Whether to consider longitudinal components of the initial vector particle explicitly or not 
    // False by default
    inline void allow_longitudinal(bool x = false)
    { 
        _longitudinal = x; 
        _helicities   = get_helicities(!x);
    };
    inline array<int,2> helicities(int n){ return _helicities[n]; };
    inline int n_helicities(){ return _helicities.size(); };

    // Momenta in different CM frames
    
    // First the ab rest frame
    inline double initial_momentum_ab(double s, double sab)
    {
        return sqrt(Kallen(sab, s, _mc2)) / sqrt(4. * sab);
    };
    inline double final_momentum_ab(double s, double sab)
    {
        return sqrt(Kallen(sab, _ma2, _mb2)) / sqrt(4. * sab);
    };

    // Next the bc rest frame
    inline double initial_momentum_bc(double s, double sbc)
    {
        return sqrt(Kallen(sbc, s, _ma2)) / sqrt(4. * sbc);
    };
    inline double final_momentum_bc(double s, double sbc)
    {
        return sqrt(Kallen(sbc, _mb2, _mc2)) / sqrt(4. * sbc);
    };

    // Finally the ac rest frame
    inline double initial_momentum_ac(double s, double sac)
    {
        return sqrt(Kallen(sac, s, _mb2)) / sqrt(4. * sac);
    };
    inline double final_momentum_ac(double s, double sac)
    {
        return sqrt(Kallen(sac, _ma2, _mc2)) / sqrt(4. * sac);
    };

    // Boundaries of  the Dalitz plot
    inline double sbc_from_sab(double s, double sab, double cos)
    {
        double thm = _ma + _mb;
        double thp = sqrt(s) - _mc;

        if (sqrt(sab) > thp || sqrt(sab) < thm) 
        {
            warning("sbc_from_sab", "Evaluating outside of physical region!");
        }

        if (abs(sqrt(sab) - thm) < EPS || abs( sqrt(sab) - thp) < EPS) return _mb2 + _mc2 - (sab - s + _mc2)*(sab - _ma2 + _mb2) / (2. * sab);

        return _mb2 + _mc2 - (sab - s + _mc2)*(sab - _ma2 + _mb2) / (2. * sab) + 2.*initial_momentum_ab(s, sab)*final_momentum_ab(s, sab)*cos;
    };

    inline double sbc_from_sac(double s, double sac, double cos)
    {
        double thm = _ma + _mc;
        double thp = sqrt(s) - _mb;

        if (sqrt(sac) > thp || sqrt(sac) < thm) 
        {
            warning("sbc_from_sac", "Evaluating outside of physical region!");
        }

        if (abs(sqrt(sac) - thm) < EPS || abs( sqrt(sac) - thp) < EPS) return _mb2 + _mc2 - (sac - s + _mb2)*(sac - _ma2 + _mc2) / (2. * sac);

        return _mb2 + _mc2 - (sac - s + _mb2)*(sac - _ma2 + _mc2) / (2. * sac) + 2.*initial_momentum_ac(s, sac)*final_momentum_ac(s, sac)*cos;
    };
    
    inline double sab_from_sbc(double s, double sbc, double cos)
    {
        double thm = _mb + _mc;
        double thp = sqrt(s) - _ma;

        if (sqrt(sbc) > thp || sqrt(sbc) < thm) 
        {
            warning("sab_from_sbc", "Evaluating outside of physical region!");
        }

        if (abs(sqrt(sbc) - thm) < EPS || abs( sqrt(sbc) - thp) < EPS) return _ma2 + _mb2 - (sbc - s + _ma2)*(sbc - _mc2 + _mb2) / (2. * sbc);
 
        return _ma2 + _mb2 - (sbc - s + _ma2)*(sbc - _mc2 + _mb2) / (2. * sbc) + 2.*initial_momentum_bc(s, sbc)*final_momentum_bc(s, sbc)*cos;
    };

    // Check whether we are considering a photon or massive vector in initial state
    bool is_photon(){ return !_longitudinal; };

    inline double Kibble(double s, double sab, double sbc)
    {
        double sac = _ma2 + _mb2 + _mc2 + s - sab - sbc;

        return sab*sbc*sac - sab*(s - _ma2)*(_mc2 - _mb2) - sbc*(s - _mb2)*(_ma2 - _mc2) 
                                          - (s*_mc2 - _mb2*_ma2)*(s -_mb2 + _ma2 - _mc2);
    };

    inline bool in_physical_region(double s, double sab, double sbc)
    {
        double phi = Kibble(s, sab, sbc);
        if (phi >= 0) return true;
        else          return false; 
    };

    // -----------------------------------------------------------------------
    private:
    
    // Whether to consider longitudinal polarizations of the initial particle explicitly
    bool _longitudinal = false;
    vector<array<int,2>> _helicities = get_helicities(true);

    // Particles can carry a string label
    string _a = "a", _b = "b", _c = "c";

    // Masses of the particles
    double _ma,  _mb,  _mc;
    double _ma2, _mb2, _mc2;
};

#endif