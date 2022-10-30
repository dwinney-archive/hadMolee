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

#include "hadronic_molecule.hpp"

namespace hadMolee
{
    // Forward declaration so we can rename ptr to kinematics as just kinematics
    // WE do this because we basically never want to work with a raw instance, but pass around a pointer
    class charmoniumlike;
    using lineshape = std::shared_ptr<hadMolee::charmoniumlike>;

    // -----------------------------------------------------------------------
    // Charmonium-like generally described a vector state that can couple to
    // a the initial-state photon. 

    // It's derived from the hadronic molecule class because we wish primarily to
    // describe the Y-mesons. However this may also be used to describe conventional
    // charmonia by setting hadronic_molecule::_molecular_coupling = 0

    class charmoniumlike : public hadronic_molecule
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor with only parameters & id (i.e. no molecular content)
        charmoniumlike( int npars, std::string id = "charmonium")
        : hadronic_molecule(), _npars(npars), _id(id)
        {};

        charmoniumlike(std::array<double,2> m, int npars, std::string id = "molecule")
        : hadronic_molecule(m), _npars(npars), _id(id)
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

    // Shortcut function to quickly create a kinematics object using smart pointers
    template<class A>
    inline std::shared_ptr<charmoniumlike> make_lineshape(std::string id)
    {
        auto model = std::make_shared<A>(id);
        return std::static_pointer_cast<charmoniumlike>(model);
    };

    // ---------------------------------------------------------------------------
    // Implementation of D1 D molecule for the Y(4260)

    class D1D_molecule : public charmoniumlike
    {
        // -----------------------------------------------------------------------
        public:

        D1D_molecule(std::string id = "Y(4260)")
        : charmoniumlike({M_D1, M_D}, 4, id)
        {
            // Set up the derivator 
            wsigma = ROOT::Math::Functor1D(this, &D1D_molecule::resigma);
            dsigma.SetFunction(wsigma);

            // Default parameters use the QQ2016 values
            set_parameters({MY_QQ2016, YBARE_QQ2016, YNM_WIDTH_QQ2016, F_Y_QQ2016});
        };

        // The propagator gains contributions from the self-energy
        std::complex<double> propagator(double s)
        {
            double E = sqrt(s);
            std::complex<double> D = E - _pole_mass - _Z*self_energy(E) + XI*_nonmol_width/2.;
            return XI * _Z / (2. * D);
        };

        // Since Y-meson is also charmonium-like it requires a photon coupling
        std::complex<double> photon_coupling()
        {
            return XI * E * _pole_mass*_pole_mass / _fY;
        };

        // When we change the pole mass, we must recalculate renormalization quantities
        // Set pole mass ,coupling, and non-mol width in a single call
        inline void set_parameters(std::vector<double> pars)
        {
            check_size(pars);

            _pole_mass              = pars[0];
            _molecular_coupling     = pars[1];
            _nonmol_width           = pars[2];
            _fY                     = pars[3];

            // Precalculate relevant quantities
            _reS  = resigma(_pole_mass);
            _redS = dsigma.Eval(_pole_mass);
            _Z    = 1. / (1. - _redS);
        };

        private:

        // Self-energy from bubble of D1 D scattering and dressed with elastic scattering
        // renomalized
        std::complex<double> self_energy(double E)
        {
            return sigma(E) - _reS - (E - _pole_mass) * _redS;
        };

        // Bare self-energy just from the bubble of D1 D scattering
        inline std::complex<double> sigma(double E)
        {
            double eps = mass_difference(E);
            double mu  = reduced_mass();
            double y   = _molecular_coupling;

            return (-XI/ (8.*PI)) * sqrt( 2.*mu*mu*mu* (eps + XI*W_D1/2.) ) * (y*y);
        };
        inline double resigma(double E){ return real(sigma(E)); };

        // Need to be able to calculate the derivative of the above self-energy
        ROOT::Math::Functor1D wsigma;
        ROOT::Math::Derivator dsigma;

        double _Z, _reS, _redS, _fY;
    };
};

#endif