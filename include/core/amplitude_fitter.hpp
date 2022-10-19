// Class which takes in an amplitude and some data and produces a fit based on chi2 minimization with MINUIT
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef FITTER
#define FITTER

#include "reaction_kinematics.hpp"
#include "amplitude.hpp"

#include <sstream> 
#include <chrono> 

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class amplitude_fitter 
{

    // --------------------------------------------------------------------
    public:

    // Constructor requires a valid amplitude to be fit
    amplitude_fitter(amplitude * amp)
    : _amplitude(amp)
    {
        // populate parameters vector of appropriate size
        for (int i = 0; i < amp->get_nParams(); i++)
        {
            _pars.push_back(i);
        };

        _minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
    };

    // Constructor requires a reaction_kinematics object
    // Optional explicit choice of minimization strategy passes to minuit object
    amplitude_fitter(amplitude * amp, string strategy, double tolerance = 1.E-6)
    : _amplitude(amp), _tolerance(tolerance)
    {
        // populate parameters vector of appropriate size
        for (int i = 0; i < amp->get_nParams(); i++)
        {
            _pars.push_back(i);
        };

        _minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy);
    };

    // Add vectors corresponding to partial-width data in a specific subchannel
    // Arguments require:
    // - specify subchannel according to labels specified in amp->_reaction_kinematics
    // - fixed sqrt(s) of the photon 
    // - vector of sqrt(sigma) where sigma is invairant mass of subsystem
    // - vector data points
    // - vector of errors
    // String id optional parameter for feeding fit results to a plotter object
    inline void add_subchannel_data(subchannel abc, double sqs, vector<double> sqsig, vector<double> data, vector<double> errors, string id = "")
    {

        // Number of data points for this set
        int n = sqsig.size();

        if (data.size() != n || errors.size() != n) 
        {
            warning("amplitude_fitter::add_subchannel_data", "Vectors received not the correct size! Skipping data set " + id + "...");
            return;
        };

        if ( id == "" ) id = subchannel_label(abc) + "_data[" + std::to_string(_subchannel_data.size()) + "]";
        data_set new_data(abc, sqs, sqsig, data, errors, id);

        _subchannel_data.push_back(new_data);

        // Add number of points to the running totals
        _N += n;
    };

    // Actually do the fit.
    void do_fit();

    // --------------------------------------------------------------------
    private:

    // Amplitude being fit
    amplitude * _amplitude;

    // MINUIT error code
    int _printLevel = 0;
    int _maxCalls   = 1E6;
    double _tolerance = 1.E-6;
    ROOT::Math::Minimizer * _minuit;
    ROOT::Math::Functor fcn;

    // minimization function
    inline double chi2(const double *pars){ return chi2_subchannels(pars); };

    // Chi2 contribution from just saved sub-channel data
    double chi2_subchannels(const double *par);

    // --------------------------------------------------------------------
    // Parameter handling

    // Container struct carrying all relevant info for each parameter involved in the fit
    struct parameter
    {
        parameter(int i)
        : _i(i), _label("par["+to_string(i)+"]")
        {};

        int _i;
        std::string _label;

        bool _fixed = false;
        double _fixed_value;

        bool _custom_limits = false;
        double _upper_limit;
        double _lower_limit;
        double _step_size = 0.1;
    };
    // Vector of size nparams holding all parameters
    std::vector<parameter> _pars;

    // convert double[] to vector<double>
    std::vector<double> convert(const double * par)
    {
        std::vector<double> result;
        for (int n = 0; n < _pars.size(); n++)
        {
            result.push_back(par[n]);
        };
        return result;
    };

    // --------------------------------------------------------------------
    // Data handling

    // Running total number of data points saved
    int _N = 0;

    // Simple container struct to hold all the relevent info for each user-added data set to fit against
    struct data_set
    {
        data_set(subchannel abc, double sqrts, vector<double> sqrtsigmas, vector<double> data, vector<double> errors, string id)
        : _subchannel(abc), _sqrts(sqrts), _sqrtsigmas(sqrtsigmas), _data(data), _errors(errors), _id(id)
        {};

        string _id;
        subchannel _subchannel;
        double _sqrts;
        vector<double> _sqrtsigmas, _data, _errors;
    };
    // Vector containing all our different data sets
    vector<data_set> _subchannel_data;

    // --------------------------------------------------------------------
    // Status message handling

    // Methods for printing out messages
    inline void divider()
    {
        cout << "--------------------------------------------------------------" << endl;
    };
    inline void new_line()
    {
        cout << endl;
    };

    // Simply return string for default particle labels
    inline string subchannel_label(subchannel x)
    {
        switch(x)
        {
            case ab: return "ab";
            case bc: return "bc";
            case ac: return "ac";
        }

        return "ERR_ID";
    };

    // Print out a summary of saved data
    void data_info();

    // Print out table of current variable status
    void variable_info(std::vector<double> pars, bool opt = 1);

    // Initialize the minimizer by feeding in saved parameter options
    void set_up(std::vector<double> starting_guess);

    // Print out summary of results
    void print_results();
};

#endif