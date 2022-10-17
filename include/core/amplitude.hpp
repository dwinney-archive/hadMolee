// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE
#define AMPLITUDE

#include "reaction_kinematics.hpp"
#include "charmoniumlike.hpp"

#include <vector>

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

// Abstract class for a generic amplitude of the e+e- -> abc process
class amplitude 
{
    // -----------------------------------------------------------------------
    public: 

    friend class amplitude_sum;

    amplitude(reaction_kinematics * xkinem, charmoniumlike * V, int n, string classname, string id = "")
    : _kinematics(xkinem), _V(V),
      _nparams(n), _id(id), _classname(classname)
    {
        array<double, 3> m = xkinem->get_masses();
        _ma  = m[0];      _mb  = m[1];      _mc  = m[2];
        _ma2 = m[0]*m[0]; _mb2 = m[1]*m[1]; _mc2 = m[2]*m[2];
    };

    // This pointer holds all kinematic information
    reaction_kinematics * _kinematics; 

    // Photon oscillates into spin-1 meson described by a charmoniumlike object
    charmoniumlike * _V;

    // Access the id tag of this amplitude
    inline void   set_id(string id){ _id = id; };
    inline string get_id()
    {
        if (_id == "") return this->_classname; 
        return this->_id;
    };

    // Debugging variable 
    inline void set_debug(int x){ _debug = x; };

    // Every amplitude needs to be able to specify how set free parameters
    virtual void set_parameters(vector<double> pars){ return; };

    // Function to output the number of parameters this amplitude takes
    inline int get_nParams(){ return _nparams; };

    // Output the amplitude at fixed helicities for photon and particle a 
    // Needs array of two helicities (inital vector and outgoing vector) 
    virtual complex<double> reduced_amplitude(int i, int j) = 0;

    // Spin-summed amplitude squared
    double probability_distribution(double s, double sab, double sbc);

    // Doubly differential partial-width
    double d2Gamma(double s, double sab, double sbc);

    // Integrated widths into given subsystem
    double dGamma_ab(double s, double sab);
    double dGamma_bc(double s, double sbc);
    double dGamma_ac(double s, double sac);

    // Fully integrated decay width
    double Gamma(double s);

    // Normalize widths to a constant N
    // This is phase-space dependent so also must specify an s 
    inline void normalize(double N, double s)
    {
        // If previously normalized, reset
        if (_normalize) 
        { 
            _normalization = 1.;
            _normalize = false;
        };

        _normalization = N / Gamma(s);
        _normalize = true;
    }
    
    // -----------------------------------------------------------------------
    protected:

    // Variable for setting cases for debugging messages
    int _debug; 

    // Quantites related to paramaters accepted
    int _nparams = 0;       // Number of parameters to expect
    vector<double> _params; // Saved parameter values
    
    // Whether or not to normalize the widths to a given constant
    bool _normalize = false;
    double _normalization = 1.;

    // Set function but kept internal (really only used by the amplitude_sum class)
    inline void set_nParams(int n){ _nparams = n; };

    // Check that passes vector is correct size for expected number of parameters
    inline void check_nParams(vector<double> params)
    {
        if (params.size() != _nparams) warning(get_id(), "Invalid number of parameters passed!");
    };

    string _id;                      // String identifier for this amplitude
    string _classname = "amplitude"; // Name of the class

    // Kinematic quantites below are saved so they do not need to be passed around inside amplitudes
    inline void update(double s, double sab, double sbc)
    {
        _W = sqrt(s); _s = s;
        _sab = sab; _sbc = sbc;
        _sac = _ma2 + _mb2 + _mc2 + s - sab - sbc;

        _updated = true;
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
    vector< vector< complex<double> > > _cached_amplitudes;
    double _cache_tolerance = EPS;
    double _cached_s, _cached_sab, _cached_sbc;
    void check_cache();

    // Total invairant energies
    double _W, _s;

    // Sub-channel energies
    double _sab, _sbc, _sac;

    // Masses of the particles
    double _ma,  _mb,  _mc;
    double _ma2, _mb2, _mc2;

    // Short cuts for characteristic angular behavior
    inline int s_wave(int i, int j) { return (i == j); }                    // S-wave is just a delta-function
    inline int d_wave(int i, int j) { return 3*(i==3)*(j==3) - (i == j); }; // D-wave 
};

// Simply amplitude with no energy dependence. 
class phase_space : public amplitude
{
    // -----------------------------------------------------------------------

    public: 

    phase_space(reaction_kinematics * xkinem, string id = "")
    : amplitude(xkinem, new charmoniumlike(""), 0, "flat_amplitude", id)
    {};

    ~phase_space()
    {
        delete _V;
    }

    // Return a constant for all the amplitudes. 
    // Since we always sum over helicities of a, we divide normalize so the spin-summed amplitude square equals 1
    // The averaging factor for the Y meson is handled in the width definition
    inline complex<double> reduced_amplitude(int i, int j) 
    { return s_wave(i,j) * XR; };
};

#endif 