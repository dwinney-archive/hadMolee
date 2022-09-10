// Abstract class for an amplitude of the e+e- -> abc reaction
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE
#define AMPLITUDE

#include "reaction_kinematics.hpp"
#include <vector>

// Abstract class for a generic amplitude of the e+e- -> abc process
class amplitude 
{
    // -----------------------------------------------------------------------
    public: 

    amplitude(reaction_kinematics * xkinem, int n, string classname, string id = "")
    : _kinematics(xkinem), _nparams(n), _id(id),
      _classname(classname)
    {};

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
    
    // -----------------------------------------------------------------------
    protected:

    // Variable for setting cases for debugging messages
    int _debug; 

    // Quantites related to paramaters accepted
    int _nparams = 0;           // Number of parameters to expect
    vector<double> _params; // Saved parameter values
    
    // Check that passes vector is correct size for expected number of parameters
    inline void check_nParams(vector<double> params)
    {
        if (params.size() != _nparams) warning(get_id(), "Invalid number of parameters passed!");
    };

    // String identifier for this amplitude
    string _id;
    string _classname = "amplitude"; // Name of the class and 
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