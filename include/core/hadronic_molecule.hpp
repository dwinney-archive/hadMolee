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

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Generic class, specific implementations for the Y and Z mesons are given below

    class hadronic_molecule
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor requires setting the masses of constituent
        // Set the constituent channel masses
        hadronic_molecule(double m1, double m2)
        : _m1(m1), _m2(m2)
        {};

        // Or given as an array
        hadronic_molecule(std::array<double,2> m)
        : _m1(m[0]), _m2(m[1])
        {};

        // Or for the "null" molecule (i.e. not a molecule) use the empty constructor
        hadronic_molecule()
        : _m1(0), _m2(0)
        {};

        // Evaluate the propagator
        // For standardization the input argument is always asssumed to be s, take the square root internally if we need E
        virtual std::complex<double> propagator(double s){ return 1.; };

        // Output the saved coupling to the constituent channel
        virtual inline double molecular_coupling(){ return _molecular_coupling; };

        // -----------------------------------------------------------------------
        protected:

        // Masses of the molecule
        double _pole_mass;

        // Width coming from decays other than m1 m2 final state
        double _nonmol_width;

        // Constituents masses 
        double _m1, _m2;
        double _molecular_coupling = 0.;

        inline double reduced_mass(){ return _m1 * _m2 / (_m1 + _m2); };
        inline double mass_difference(double E){ return E - _m1 - _m2; };

        // Intermediate state threshold openings of the constituent channel
        inline double Wth(){ return _m1 + _m2; };
        inline double sth(){ return Wth()*Wth(); };

        // Self-energy loop function
        virtual std::complex<double> self_energy(double x){ return 0.; };
    };

    // ---------------------------------------------------------------------------
    // D* D molecule relevant for the Z meson

    class DsD_molecule : public hadronic_molecule
    {
        // -----------------------------------------------------------------------
        public:

        DsD_molecule()
        : hadronic_molecule(M_DSTAR, M_D)
        {
            // Mass and Width from PDG
            _pole_mass          = M_ZC3900;
            _total_width        = W_ZC3900;
            
            // Coupling taken from [1]
            _molecular_coupling      = ZBARE_QQ2016;  
            
            // residual width taken to recover the full PDG width at the pole
            _nonmol_width = _total_width - 2.* imag(self_energy(_pole_mass));
        };

        // The propagator gains contributions from the self-energy
        inline std::complex<double> propagator(double s)
        {
            double E = sqrt(s);
            double z = _molecular_coupling;

            std::complex<double> D = E - _pole_mass + XI * (z*z*self_energy(E) + _nonmol_width/2.);
            
            return XI / (2.*D);
        };  

        // Self-energy from bubble diagram of D* D scattering 
        inline std::complex<double> self_energy(double E)
        {
            double eps = mass_difference(E);
            double mu  = reduced_mass();
            
            return (1. / (8.*PI)) * sqrt(2.*mu*mu*mu*std::abs(eps)) * ( XR*(eps>=0) + XI*(eps<0) );
        };

        // -----------------------------------------------------------------------
        private:

        // Total width of the Z from the PDG
        double _total_width;
    };
};

#endif