// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE
#define AMPLITUDE

#include "reaction_kinematics.hpp"

#include <vector>

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

// Abstract class for a generic amplitude of the e+e- -> abc process
class amplitude 
{
    // -----------------------------------------------------------------------
    public: 

    amplitude(reaction_kinematics * xkinem, int n, string classname, string id = "")
    : _kinematics(xkinem), _nparams(n), _id(id),
      _classname(classname)
    {
        array<double, 3> m = xkinem->get_masses();
        _ma  = m[0];      _mb  = m[1];      _mc  = m[2];
        _ma2 = m[0]*m[0]; _mb2 = m[1]*m[1]; _mc2 = m[2]*m[2];
    };

    // This pointer holds all kinematic information
    reaction_kinematics * _kinematics; 

    // Access the id tag of this amplitude
    inline void   set_id(string id){ _id = id; };
    inline string get_id()
    {
        if (_id == "") return this->_classname; 
        return this->_id;
    };

    // Debugging variable 
    inline void set_debug(int x){ _debug = x; };

    // Output the amplitude at fixed helicities for photon and particle a 
    // Needs array of two helicities (inital vector and outgoing vector) 
    // The total invariant mass, s (mass of decaying particle)
    // and we choose sigma_ab and sigma_bc to be the independent variables
    virtual complex<double> helicity_amplitude(array<int,2> helicities, double s, double sab, double sbc) = 0;
    inline  complex<double> helicity_amplitude(int n, double s, double sab, double sbc)
    {
        return helicity_amplitude( get_helicities(n, _kinematics->is_photon()), s, sab, sbc);
    };


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
    
    // -----------------------------------------------------------------------
    protected:

    // Variable for setting cases for debugging messages
    int _debug; 

    // Quantites related to paramaters accepted
    int _nparams = 0;       // Number of parameters to expect
    vector<double> _params; // Saved parameter values
    
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
    };

    // Total invairant energies
    double _W, _s;

    // Sub-channel energies
    double _sab, _sbc, _sac;

    // Masses of the particles
    double _ma,  _mb,  _mc;
    double _ma2, _mb2, _mc2;
};

// Simply amplitude with no energy dependence. 
class flat_amplitude : public amplitude
{
    // -----------------------------------------------------------------------
    public: 
    flat_amplitude(reaction_kinematics * xkinem, string id = "")
    : amplitude(xkinem, 0, "flat_amplitude", id)
    {};

    // Return a constant for all the amplitudes. 
    // Since we always sum over helicities of a, we divide normalize so the spin-summed amplitude square equals 1
    // The averaging factor for the Y meson is handled in the width definition
    inline complex<double> helicity_amplitude(array<int,2> helicities, double s, double sab, double sbc) 
    { return XR / sqrt(3.); };

    // -----------------------------------------------------------------------
    protected:
};

#endif 