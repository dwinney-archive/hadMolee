// Container class holding all quantities relevant for a hadronic molecule state
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// References:

// [1] arXiv:1310.2190 [hep-ph]
// ---------------------------------------------------------------------------


#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include "Math/Functor.h"
#include "Math/Derivator.h"
#include "constants.hpp"

namespace hadMolee
{
    class molecular;

    // The molecule pieces are not assumed to be the degrees of freedom of a given model by themselves
    // These are either contained in the lineshape or amplitude classes and do not need to be shared
    using molecule = std::shared_ptr<molecular>;
    
    // Shortcut functions to quickly create a kinematics object using smart pointers

    // In case the constructor requires constituent particle masses
    inline molecule make_molecule(double a, double b)
    {
        molecule model = std::make_shared<molecular>(a, b);
        return model;
    };

    // But in most cases we just have those preset in a model and require no arguments
    template<class A>
    inline molecule make_molecule()
    {
        molecule model = std::make_shared<A>();
        return model;
    };

    // If we have a lineshape model that also contains a molecular piece
    // this method extracts it
    template<class A>
    inline molecule get_molecular_component(A m)
    {
        molecule ptr = std::dynamic_pointer_cast<molecular>(m);
        return ptr;
    };

    // This simply checks if a given model has a molecular component or not
    template<class A>
    inline bool is_molecular(A m)
    {
        molecule ptr = std::dynamic_pointer_cast<molecular>(m);
        if (ptr) return true;
        return false;
    };

    // ---------------------------------------------------------------------------
    // Generic class, specific implementations for the Y and Z mesons are given below

    class molecular
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor requires setting the masses of constituent
        // Set the constituent channel masses
        molecular(double m1, double m2)
        : _m1(m1), _m2(m2)
        {};

        // Or given as an array
        molecular(std::array<double,2> m)
        : _m1(m[0]), _m2(m[1])
        {};

        // Or given as an array
        molecular(std::array<double,2> m, std::array<double,2> w)
        : _m1(m[0]), _m2(m[1]), _w1(w[0]), _w2(w[1])
        {};

        virtual void set_parameters(std::array<double,3> pars)
        {
            _mass     = pars[0];
            _nm_width = pars[1];
            _coupling = pars[2];

            _reSigmaPole = std::real(self_energy(_mass*_mass));
        };

        // Evaluate the propagator
        // For standardization the input argument is always asssumed to be s, take the square root internally if we need E
        virtual complex propagator(double s)
        { 
            return 1./(sqrt(s)-_mass-_coupling*_coupling*(self_energy(s)-_reSigmaPole)+I*_nm_width/2.)/(2*_mass);
        };

        // Output the saved coupling to the constituent channel
        virtual inline double coupling(){ return _coupling; }; 
        virtual inline double mass()    { return _mass;     };  

        // Self-energy loop function from DR
        // Always a function of s (GeV^2)
        virtual inline complex self_energy(double s)
        {
            double mu   = .770;  // renormalization at rho mass
            double amu  = 0.;  // DR renomalization coefficient

            complex m1  = _m1 + I*_w1/2.;
            complex m2  = _m2 + I*_w2/2.;

            complex rho = csqrt(Kallen(s+I*0, m1*m1, m2*m2))/s;
            complex xi  = 1. - (m1+m2)*(m1+m2)/s;

            // Imaginary part comes from this piece
            complex logs = rho*log((xi + rho)/(xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1);
            
            // To match the energy dependence of the dimensionally regularization formula,
            // we add this piece
            complex  DR = amu + log(m1*m1/mu*mu) + (s-m1*m1+m2*m2)/s*log(m2/m1);

            return (DR + logs)/(16.*PI*PI);
        };

        // Return masses 
        inline std::array<double,2> constituent_masses(){ return {_m1, _m2}; };

        double _reSigmaPole = 0;
        // -----------------------------------------------------------------------
        protected:

        // Masses of the molecule
        double _mass;

        // Width coming from decays other than m1 m2 final state
        double _nm_width = 0;

        // Constituents masses 
        double _m1 = 0, _m2 = 0;
        double _w1 = 0, _w2 = 0;
        double _coupling = 0.;

        inline double reduced_mass(){ return _m1 * _m2 / (_m1 + _m2); };
        inline double mass_difference(double E){ return E - _m1 - _m2; };

        // Intermediate state threshold openings of the constituent channel
        inline double Wth(){ return _m1 + _m2; };
        inline double sth(){ return Wth()*Wth(); };
    };
};

#endif