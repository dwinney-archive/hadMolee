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
#include <memory>

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

namespace hadMolee
{
    // Forward declaration so we can rename ptr to kinematics as just kinematics
    // WE do this because we basically never want to work with a raw instance, but pass around a pointer
    class reaction_kinematics;
    using kinematics = std::shared_ptr<hadMolee::reaction_kinematics>;

    // Generic particle labels
    enum particle{a, b, c};

    // Three possible subchannels following the labelings inside amp->_reaction_kinematics
    enum subchannel{ab, ba = ab, bc, cb = bc, ac, ca = ac};

    // Abstract class for a generic amplitude of the e+e- -> abc process
    class reaction_kinematics 
    {
        // -----------------------------------------------------------------------
        public: 

        // Constructor with only masses and number of parameters
        reaction_kinematics(std::array<double,3> m)
        {
            set_particle_masses(m);
        };

        // Constructor that also allows labels
        reaction_kinematics(std::array<double,3> m, std::array<std::string,3> labels)
        {
            set_particle_masses(m);
            set_particle_labels(labels);
        };

        // Return the string identifier
        inline std::string get_id(){ return "[" + particle_a() + ", " + particle_b() + ", " + particle_c() +  "]"; };

        // -----------------------------------------------------------------------
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

            if (std::abs(sqrt(sab) - thm) < EPS || std::abs( sqrt(sab) - thp) < EPS) return _mb2 + _mc2 - (sab - s + _mc2)*(sab - _ma2 + _mb2) / (2. * sab);

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

            if (std::abs(sqrt(sac) - thm) < EPS || std::abs( sqrt(sac) - thp) < EPS) return _mb2 + _mc2 - (sac - s + _mb2)*(sac - _ma2 + _mc2) / (2. * sac);

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

        inline double sac_from_sab(double s, double sab, double cos)
        {
            double sbc = sbc_from_sab(s, sab, cos);
            return _ma2 + _mb2 + _mc2 + s - sbc - sab;
        };  

        inline double sac_from_sbc(double s, double sbc, double cos)
        {
            double sab = sab_from_sbc(s, sbc, cos);
            return _ma2 + _mb2 + _mc2 + s - sab - sbc;
        };

        inline double Kibble(double s, double sab, double sbc)
        {
            double sac = _ma2 + _mb2 + _mc2 + s - sab - sbc;
            warning("Kibble was weird last time I checked use with caution");
            return sab*sbc*sac - sab*(s - _ma2)*(_mc2 - _mb2) - sbc*(s - _mb2)*(_ma2 - _mc2) 
                                            - (s*_mc2 - _mb2*_ma2)*(s -_mb2 + _ma2 - _mc2);
        };

        inline bool in_physical_region(double s, double sab, double sbc)
        {
            return (sbc >= sbc_from_sab(s, sab, -1.)) && (sbc <= sbc_from_sab(s, sab, +1.));
        };
        inline bool in_physical_region(double s, double sig, subchannel x)
        {
            double upper, lower;
            switch (x)
            {
                case ab: 
                {
                    upper = (sqrt(s) - _mc);
                    lower = (_ma + _mb);
                    break;
                }
                case bc: 
                {
                    upper = (sqrt(s) - _ma);
                    lower = (_mb + _mc);
                    break;
                }
                case ac: 
                {
                    upper = (sqrt(s) - _mb);
                    lower = (_ma + _mc);
                    break;
                }
                default: return false;
            };

            return ((sqrt(sig) <= upper) && (sqrt(sig) >= lower));
        };

        // Set particle labels (useful for plotting)
        inline void set_particle_labels(std::string a, std::string b, std::string c){ _a = a; _b = b; _c = c; };
        inline void set_particle_labels(std::array<std::string,3> l){ set_particle_labels(l[0], l[1], l[2]); };
        inline std::string particle_a(){ return _a; };
        inline std::string particle_b(){ return _b; };
        inline std::string particle_c(){ return _c; };

        inline std::string subchannel_label(subchannel x)
        {
            switch (x) 
            {
                case ab: return particle_a() + " " + particle_b();
                case ac: return particle_a() + " " + particle_c();
                case bc: return particle_b() + " " + particle_c();
                default:
                {
                    warning("reaction_kinematics::subchannel_label()", "Subchannel not found! Returning empty string...");
                    return "";
                }
            };

            return "";
        };

        inline void set_particle_masses(double ma, double mb, double mc)
        { 
            _ma  = ma;     _mb = mb;    _mc  = mc;
            _ma2 = ma*ma; _mb2 = mb*mb; _mc2 = mc*mc;
        };
        inline void set_particle_masses(std::array<double,3> m)
        {
            set_particle_masses(m[0], m[1], m[2]);
        };
        inline std::array<double,3> get_masses(){ return {_ma, _mb, _mc}; };
        inline double mass_a(){ return _ma; };
        inline double mass_b(){ return _mb; };
        inline double mass_c(){ return _mc; };

        // -----------------------------------------------------------------------
        // Momenta in the lab frame

        // Particle a momentum in the vector helicity frame
        inline double decay_momentum_a(double s, double sbc)
        {
            return sqrt(Kallen(s, sbc, _ma2)) / (2. * sqrt(s));
        };

        // Particle b momentum in the vector helicity frame
        inline double decay_momentum_b(double s, double sac)
        {
            return sqrt(Kallen(s, sac, _mb2)) / (2. * sqrt(s));
        };

        // Particle c momentum in the vector helicity frame
        inline double decay_momentum_c(double s, double sab)
        {
            return sqrt(Kallen(s, sab, _mc2)) / (2. * sqrt(s));
        };
        
        // -----------------------------------------------------------------------
        // Relative angles
        
        // Polar angle between b and c in the lab frame
        inline double cos_bc(double s, double sab, double sbc, double sac)
        {
            double pb = decay_momentum_b(s, sac);
            double pc = decay_momentum_c(s, sab);

            double Eb = sqrt(pb*pb + _mb2);
            double Ec = sqrt(pc*pc + _mc2);

            return (2.*Eb*Ec - sbc + _mb2 + _mc2) / (2.*pb*pc);
        };

        // -----------------------------------------------------------------------
        
        protected:

        // String identifier
        std::string _id;

        // Variable for setting cases for debugging messages
        int _debug; 

        // Masses of the particles
        double _ma,  _mb,  _mc;
        double _ma2, _mb2, _mc2;

        // Particles can carry a string label
        std::string _a = "a", _b = "b", _c = "c";
    };

    // Shortcut function to quickly create a kinematics object using smart pointers
    inline kinematics make_kinematics(std::array<double,3> m)
    {
        return std::make_shared<reaction_kinematics>(m);
    };

    // Shortcut function to quickly create a kinematics object using smart pointers
    inline kinematics make_kinematics(std::array<double,3> m, std::array<std::string,3> labels)
    {
        return std::make_shared<reaction_kinematics>(m, labels);
    };
};

#endif 