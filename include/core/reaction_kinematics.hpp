// Abstract class for cross-section describing the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef KINEMATICS
#define KINEMATICS

#include "constants.hpp"

#include <array>
#include <vector>

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

// Abstract class for a generic amplitude of the e+e- -> abc process
class reaction_kinematics 
{
    // -----------------------------------------------------------------------
    public: 

    // Constructor with only masses and number of parameters
    reaction_kinematics(array<double,3> m)
    {
        set_particle_masses(m);
    };

    // Constructor that also allows labels
    reaction_kinematics(array<double,3> m, array<string,3> labels)
    {
        set_particle_masses(m);
        set_particle_labels(labels);
    };

    // ------------------------------------------------------------------------------------------------------------------
    // Kinematic quantities 

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

    // Set particle labels (useful for plotting)
    inline void set_particle_labels(string a, string b, string c){ _a = a; _b = b; _c = c; };
    inline void set_particle_labels(array<string,3> l){ set_particle_labels(l[0], l[1], l[2]); };
    inline string particle_a(){ return _a; };
    inline string particle_b(){ return _b; };
    inline string particle_c(){ return _c; };

    inline void set_particle_masses(double ma, double mb, double mc)
    { 
         _ma  = ma;     _mb = mb;    _mc  = mc;
         _ma2 = ma*ma; _mb2 = mb*mb; _mc2 = mc*mc;
    };
    inline void set_particle_masses(array<double,3> m)
    {
        set_particle_masses(m[0], m[1], m[2]);
    };
    inline array<double,3> get_masses(){ return {_ma, _mb, _mc}; };
    inline double mass_a(){ return _ma; };
    inline double mass_b(){ return _mb; };
    inline double mass_c(){ return _mc; };

    // -----------------------------------------------------------------------
    // Vectors and tensors

    // Lepton-tensor gives the e+ e- -> photon interaction
    complex<double> production_tensor(int i, int j, double s)
    {
        // In massless approximation, this reduces to the polarization sum of transverse vector meson
        complex<double> lepton_tensor = E*E * 2.*s * XR*((i == j) - (i==3)*(j==3));
        complex<double> photon_propagator = norm(XI / s);
        return lepton_tensor  * photon_propagator;
    };

    // -----------------------------------------------------------------------
    protected:

    // Variable for setting cases for debugging messages
    int _debug; 

    // Masses of the particles
    double _ma,  _mb,  _mc;
    double _ma2, _mb2, _mc2;

    // Particles can carry a string label
    string _a = "a", _b = "b", _c = "c";
};

#endif 