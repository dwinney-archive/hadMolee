// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE
#define AMPLITUDE

#include "kinematics.hpp"
#include "lineshape.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

namespace hadMolee
{
    // Forward declare amplitude for typedefs below
    class amplitude_base;

    // We really only want to work with pointers to amplitudes and
    // almost never with amplitudes themselves
    using amplitude  = std::shared_ptr<amplitude_base>;

    // Key class with private constructor, this is required to initialize an amplitude
    // but itself may only be initilaized by friend functions
    class amplitude_key 
    {
        private:
        amplitude_key(){};

        template<class A>
        friend amplitude make_amplitude(kinematics, lineshape, std::string);
        template<class A>
        friend amplitude make_amplitude(kinematics, std::string);

        friend amplitude operator+(amplitude, amplitude);
    };
    
    // ---------------------------------------------------------------------------
    // These methods deal with manipulating amplitudes in shared_ptr form

    bool are_compatible(amplitude a, amplitude b);

    // Smart pointer "constructor"
    template<class A>
    amplitude make_amplitude(kinematics xkinem, lineshape V, std::string id)
    {
        auto amp = std::make_shared<A>(amplitude_key(), xkinem, V, id);
        return std::static_pointer_cast<amplitude_base>(amp);
    };

    template<class A>
    amplitude make_amplitude(kinematics xkinem, std::string id)
    {
        auto amp = std::make_shared<A>(amplitude_key(), xkinem, id);
        return std::static_pointer_cast<amplitude_base>(amp);
    };

    // Smart pointer "constructor" for a sum by adding two amplitudes together
    amplitude operator+(amplitude a, amplitude b);

    // Internally we can alias amplitude::add() with the += operator
    void operator+=(amplitude a,  amplitude b);

    // ---------------------------------------------------------------------------
    // Abstract class for a generic amplitude of the e+e- -> abc process
    class amplitude_base
    {
        // -----------------------------------------------------------------------
        public: 

        // Public constructors however these shouldnt be used directly
        // Rather the make_amplitude<>() method should be treated like the 'constructor'
        amplitude_base(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "")
        : _kinematics(xkinem), _V(V),
        _nparams(0), _id(id), _classname("amplitude")
        {
            set_up(xkinem);
        };

        // Constructor which specified number of parameters and a default classname string
        amplitude_base(amplitude_key key, kinematics xkinem, lineshape V, int n, std::string classname, std::string id = "")
        : _kinematics(xkinem), _V(V),
        _nparams(n), _id(id), _classname(classname)
        {
            set_up(xkinem);
        };

        // This pointer holds all kinematic information
        kinematics _kinematics; 

        // Photon oscillates into spin-1 meson described by a charmoniumlike object
        lineshape _V;

        // Access the id tag of this amplitude
        inline void        set_id(std::string id){ _id = id; };
        inline std::string get_id()
        {
            if (_id == "") return this->_classname; 
            return this->_id;
        };

        // Debugging variable 
        inline void set_debug(int x){ _debug = x; };

        // Allocate an aggragated vector of parameters to individual amplitudes
        virtual void set_parameters(std::vector<double> x);

        // Function to output the number of parameters this amplitude takes
        inline int N_parameters(){ return _nparams; };

        // Output the amplitude at fixed helicities for photon and particle a 
        // Needs array of two helicities (inital vector and outgoing vector) 
        virtual complex reduced_amplitude(cartesian_index i, cartesian_index j);

        // The above needs to be supplied by the user
        // This version is used internally and just automatically provides the _update and recalcualte() check
        inline complex reduced_amplitude_checked(cartesian_index i, cartesian_index j)
        {
            // IF our kinematics has been changed, recalculate relevant quantities
            if (updated()) recalculate();
            return reduced_amplitude(i, j);
        };

        // 2->3 distribution including a polar orientation angle to the beam
        virtual double decay_distribution(double s, double sab, double sbc, double cos);
        virtual double decay_distribution(double s, double sab, double sbc);

        // Doubly differential partial-width
        double differential_xsection(double s, double sab, double sbc);

        // Integrated widths into given subsystem
        double differential_xsection(subchannel chan, double s, double sigma);

        // Integrated production cross-section
        double integrated_xsection(double S);

        // Multiply decay widths with an arbitrary constant
        inline void normalize(double N)
        {
            // If previously normalized, reset
            if (_normalize) 
            { 
                _normalization = 1.;
                _normalize = false;
            };

            _normalization = N;
            _normalize = true;
        }
        
        // -----------------------------------------------------------------------
        protected:

        // ----------------------------------------
        // These first memebers are related to being able to sum amplitudes together

        // Sums hold a vector of pointers to the constituent amplitudes
        bool is_sum(){ return (_sub_amps.size() > 0);  };
        std::vector<amplitude> _sub_amps;

        // Make sure all amplitudes being summed are compatible
        // This is efffectively the same as the are_compatible() method except comparisions are made against the interally stored instances 
        bool is_compatible(amplitude new_amp);
    
        // Add a new amplitude to the vector
        inline void add(amplitude new_amp)
        {
            if (is_compatible(new_amp))
            {
                _sub_amps.push_back(new_amp);
                set_nParams(N_parameters() + new_amp->N_parameters());
            };
        };

        // Add a whole vector of new amplitudes
        inline void add(std::vector<amplitude> new_amps)
        {
            for (auto amp : new_amps)
            {
                add(amp);
            };
        };

        // Need to make binary operators able to access the _sub_amps vector
        friend void add(amplitude new_amp);
        friend void add(amplitude new_amps);
        friend void operator+=(amplitude a, amplitude b);
        friend amplitude operator+(amplitude a, amplitude b);
    
        // -----------------------------------------------------------------------
        // Setters and debugging 

        // Variable for setting cases for debugging messages
        int _debug; 

        // Quantites related to paramaters accepted
        int _nparams = 0;       // Number of parameters to expect
        std::vector<double> _params; // Saved parameter values
        
        // Whether or not to normalize the widths to a given constant
        bool   _normalize = false;
        double _normalization = 1.;

        // Set function but kept internal (really only used by the amplitude_sum class)
        inline void set_nParams(int n){ _nparams = n; };

        // Check that passes vector is correct size for expected number of parameters
        inline void check_nParams(std::vector<double> params)
        {
            if (params.size() != _nparams) warning(get_id(), "Invalid number of parameters passed!");
        };

        std::string _id = "";                 // String identifier for this amplitude
        std::string _classname = "amplitude"; // Name of the class

        // Kinematic quantites below are saved so they do not need to be passed around inside amplitudes
        inline void set_up(kinematics xkinem)
        {
            std::array<double, 3> m = xkinem->get_masses();
            _ma  = m[0];      _mb  = m[1];      _mc  = m[2];
            _ma2 = m[0]*m[0]; _mb2 = m[1]*m[1]; _mc2 = m[2]*m[2];
        };

        // -----------------------------------------------------------------------
        // Related to caching 

        // Total invairant energies
        double _W, _s;

        // Orientation polar angle to the beam
        double _cos;

        // Sub-channel energies
        double _sab, _sbc, _sac;

        // Masses of the particles
        double _ma,  _mb,  _mc;
        double _ma2, _mb2, _mc2;

        // Caching energy variables
        inline void update(double s, double sab, double sbc, double cos)
        {
            _W   = sqrt(s); _s   = s; _cos = cos;
            _sab = sab;     _sbc = sbc;
            _sac = _ma2 + _mb2 + _mc2 + s - sab - sbc;

            _updated = true;

            // Check if we need to update our precalculated amplitudes 
            check_decay_cache();
        };

        // If recently updated, reset the _updated flag and return 1
        // if not then return 0
        // This is used to check if quantities need to be recalculated if energies have not changed
        bool _updated = false;
        inline bool updated()
        {
            if (_updated)
            {
                _updated = false;
                return true;
            }
            return false;
        }; 

        // If updated() == true, this function should be called to recalculate all member data of the amplitude
        // which depends on the kinematic state
        virtual void recalculate(){ return; };

        // Save all the components of the reduced amplitude at each step of sab and bc to avoid having to recalculate
        std::array<std::array<double,3>,3> _cached_decay_tensor;
        double _cache_tolerance = EPS;
        double _cached_s, _cached_sab, _cached_sbc, _cached_cos;
        void check_decay_cache();

        // Short cuts for characteristic angular behavior
        inline double phat(cartesian_index i)
        {
            switch (i)
            {
                case x : return sqrt(1. - _cos*_cos);
                case y : return 0;
                case z : return _cos;
            };
        };

        // Differential widths in terms of specific channels
        double dGamma_ab(double s, double sab);
        double dGamma_bc(double s, double sbc);
        double dGamma_ac(double s, double sac);
    };

    // Simply amplitude with no energy dependence. 
    class phase_space : public amplitude_base
    {
        // -----------------------------------------------------------------------

        public: 

        phase_space(amplitude_key key, kinematics xkinem, lineshape V, std::string id = "")
        : amplitude_base(key, xkinem, V, 0, "flat_amplitude", id)
        {};

        phase_space(amplitude_key key, kinematics xkinem, std::string id = "")
        : amplitude_base(key, xkinem, nullptr, 0, "flat_amplitude", id)
        {};

        // Return a constant for all the amplitudes. 
        // Since we always sum over helicities of a, we divide normalize so the spin-summed amplitude square equals 1
        // The averaging factor for the Y meson is handled in the width definition
        inline double decay_distribution(double s, double sab, double sbc, double cos)
        {
            return 1;
        }
    };
};

#endif 