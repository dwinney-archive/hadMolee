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
#include "debug.hpp"

namespace hadMolee
{
    // Forward declaration so we can rename ptr to kinematics as just kinematics
    // WE do this because we basically never want to work with a raw instance, but pass around a pointer
    class charmoniumlike;

    // A line-shape is a model for something that is charmonioumlike that can be passed around
    // i.e. a pointer to charmoniumlike object
    using lineshape = std::shared_ptr<hadMolee::charmoniumlike>;

    // Shortcut function to quickly create a kinematics object using smart pointers
    template<class A>
    inline lineshape make_lineshape(std::string id)
    {
        auto model = std::make_shared<A>(id);
        return std::static_pointer_cast<charmoniumlike>(model);
    };

    // -----------------------------------------------------------------------
    // Charmonium-like generally described a vector state that can couple to
    // a the initial-state photon. 

    // It's derived from the hadronic molecule class because we wish primarily to
    // describe the Y-mesons. However this may also be used to describe conventional
    // charmonia by setting hadronic_molecule::_molecular_coupling = 0

    class charmoniumlike
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor with only a string id 
        charmoniumlike(std::string id = "e+e-")
        : _npars(0), _id(id)
        {};

        // Constructor that also sets nonzero number of free parameters
        charmoniumlike( int npars, std::string id = "e+e-")
        : _npars(npars), _id(id)
        {};

        virtual std::complex<double> propagator(double s){ return 1.; };
        virtual std::complex<double> photon_coupling(){    return 1.; };

        // String identifier
        inline std::string get_id(){ return _id; };

        // Set pole mass ,coupling, and non-mol width in a single call
        virtual inline void set_parameters(std::vector<double> pars)
        {
            check_size(pars);
            return;
        };  

        // Access number of free parameters from outside
        inline int N_parameters(){ return _npars; };

        // -----------------------------------------------------------------------
        protected:
        
        // Name identifier
        std::string _id;

        int _npars = 0;
        void check_size(std::vector<double> pars)
        {
            if (pars.size() != _npars)
            {
                warning("lineshape", "Wrong number of parameters given! Expected 3 but recieved " + std::to_string(pars.size()) + ". Results may vary...");
            };
        }
    };
};

#endif