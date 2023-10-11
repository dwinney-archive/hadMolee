// Container class holding all quantities relevant for a charmoniumlike particle
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// References:

// [1] arXiv:1310.2190 [hep-ph]
// ---------------------------------------------------------------------------


#ifndef LINESHAPE_HPP
#define LINESHAPE_HPP

#include "constants.hpp"

namespace hadMolee
{
    // Forward declaration so we can rename ptr to lineshape
    // WE do this because we basically never want to work with a raw instance, but pass around a pointer
    class charmonium;

    // A line-shape is a model for something that is charmoniumthat can be passed around
    // i.e. a pointer to charmoniumobject
    using lineshape = std::shared_ptr<hadMolee::charmonium>;

    // Shortcut function to quickly create a  lineshape object using smart pointers
    template<class A>
    inline lineshape make_lineshape()
    {
        auto model = std::make_shared<A>();
        return std::static_pointer_cast<charmonium>(model);
    };

    // Shortcut function to quickly create a kinematics object using smart pointers
    template<class A, class B>
    inline lineshape make_lineshape(B args)
    {
        auto model = std::make_shared<A>(args);
        return std::static_pointer_cast<charmonium>(model);
    };

    // -----------------------------------------------------------------------
    // Charmonium-like generally described a vector state that can couple to
    // a the initial-state photon. 

    // It's derived from the hadronic molecule class because we wish primarily to
    // describe the Y-mesons. However this may also be used to describe conventional
    // charmonia by setting hadronic_molecule::_molecular_coupling = 0

    class charmonium
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor with only a string id 
        charmonium()
        : _npars(3)
        {};

        charmonium(std::vector<double> pars)
        : _npars(3)
        {
            set_parameters(pars);
        };

        // Constructor that also sets nonzero number of free parameters
        charmonium(int npars)
        : _npars(npars)
        {};

        // Defaults to a relativistic Breit-Wigner but can be overridden by a lineshape model
        virtual complex propagator(double s){ return 1./(s - _mass*_mass + I*_mass*_width); };

        // Photon coupling always takes the VMD form in terms of the decay constant
        inline double photon_coupling(){ return E * _mass*_mass / _fV; };
        inline double mass()           { return _mass; };

        // String identifier
        inline void        set_id(std::string x ){ _id = x;    };
        inline std::string get_id()              { return _id; };

        // Set pole mass ,coupling, and non-mol width in a single call
        virtual inline void set_parameters(std::vector<double> pars)
        {
            check_size(pars);
            _mass  = pars[0];
            _width = pars[1];
            _fV    = decay_constant(pars[2]);
            return;
        };  

        // Access number of free parameters from outside
        inline int N_parameters(){ return _npars; };

        // -----------------------------------------------------------------------
        protected:

        // Calculate the decay constant from to  V -> e+ e- branching fraction
        inline double decay_constant(double BR){ return sqrt( 4.*PI*ALPHA* ALPHA * _mass / (3.*_width*BR)); };
        
        // Name identifier
        std::string _id;

        // Saved properties include mass, width and decay constant
        double _mass = 0., _width = 0., _fV = 0.;

        int _npars = 0;
        void check_size(std::vector<double> pars)
        {
            if (pars.size() != _npars)
            {
                warning("lineshape", "Wrong number of parameters given! Expected " + std::to_string(_npars) + " but recieved " + std::to_string(pars.size()) + ". Results may vary...");
            };
        }
    };
};

#endif