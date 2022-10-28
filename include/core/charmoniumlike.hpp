// Container class holding all quantities relevant for a charmoniumlike particle
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------
// References:

// [1] arXiv:1310.2190 [hep-ph]
// ---------------------------------------------------------------------------


#ifndef CHARMONIUMLIKE
#define CHARMONIUMLIKE

#include "constants.hpp"

namespace hadMolee
{
    class charmoniumlike
    {
        // -----------------------------------------------------------------------
        public:

        // Empty constructor
        charmoniumlike( int npars, std::string id = "charmonium")
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
                warning("charmoniumlike", "Wrong number of parameters given! Expected 3 but recieved " + std::to_string(pars.size()) + ". Results may vary...");
            };
        }
    };
};

#endif