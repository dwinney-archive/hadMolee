// Class which takes in an amplitude and some data and produces a fit based on chi2 minimization with MINUIT
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef FITTER
#define FITTER

#include "reaction_kinematics.hpp"
#include "amplitude.hpp"
#include "data_formatting.hpp"

#include <sstream> 
#include <chrono> 

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace hadMolee
{
    class amplitude_fitter 
    {

        // --------------------------------------------------------------------
        public:

        // Constructor requires a valid amplitude to be fit
        amplitude_fitter(amplitude * amp)
        : _amplitude(amp), _V(amp->_V)
        {
            // Extract how many parameters we should expect
            _N_V    = amp->_V->N_parameters();  // # of parameters from production & lineshape of vector meson
            _N_amp  = amp->N_parameters();      // # of parameters from decay to specific final state
            _N_pars = _N_V + _N_amp;

            // populate parameters vector of appropriate size
            for (int i = 0; i < _N_pars; i++)
            {
                _pars.push_back(i);
            };

            _minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
        };

        // Constructor requires a reaction_kinematics object
        // Optional explicit choice of minimization strategy passes to minuit object
        amplitude_fitter(amplitude * amp, std::string strategy, double tolerance = 1.E-6)
        : _amplitude(amp), _tolerance(tolerance)
        {
            // Extract how many parameters we should expect
            int _N_V    = amp->_V->N_parameters();  // # of parameters from production & lineshape of vector meson
            int _N_amp  = amp->N_parameters();      // # of parameters from decay to specific final state
            int _N_pars = _N_V + _N_amp;

            // populate parameters vector of appropriate size
            for (int i = 0; i < _N_pars; i++)
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
        void add_subchannel_data(subchannel abc, double sqs, std::vector<double> sqsig, std::vector<double> data, std::array<std::vector<double>,2> errors, std::string id = "");

        // Alternatively you can jsut specify the subchannel, fixed center-of-mass energy, and point to a data file
        void add_subchannel_data(subchannel abc, double sqs, std::string filename, std::string id = "");

        // Give each parameter a string name for outputting results
        void set_parameter_labels(std::vector<std::string> labels);

        // Set limits of a parameter specificed by index or name 
        void set_parameter_limits(int i, std::array<double,2> ranges, double step = 0.1);
        inline void set_parameter_limits(std::string name, std::array<double,2> ranges, double step = 0.1) { set_parameter_limits(find_parameter(name), ranges, step); };

        // Specify a parameter should be considered fix in the fit loop
        inline void fix_parameter(int i) { _pars[i]._fixed = true; };
        inline void fix_parameter(std::string name){ fix_parameter(find_parameter(name)); };

        // Unfix a parameter
        inline void free_parameter(int i){ _pars[i]._fixed = false; };

        //Utility to change print level in TMinuit, default is to surpress all messages
        inline void set_print_level(int n){ _printLevel = n; };

        // Change the maximal number of function calls allowd 
        inline void set_max_calls(int n){ _maxCalls = n; };

        // Actually do the fit.
        void do_fit(std::vector<double> starting_guess);

        // --------------------------------------------------------------------
        private:

        // Decay amplitude being fit 
        // This describes the definite final state we are considering. 
        // In future this can be a vector to include multiple final states with their own data sets
        amplitude * _amplitude;

        // This pointer is the shared vector meson object that all the amplitudes being fit share
        // This can be our charmonium or Y-meson describing the total center-of-mass energy dependence
        charmoniumlike * _V;

        // MINUIT error code
        int _printLevel   = 0;
        int _maxCalls     = 1E6;
        double _tolerance = 1.E-6;
        ROOT::Math::Minimizer * _minuit;
        ROOT::Math::Functor fcn;

        // minimization function
        inline double chi2(const double *pars){ return chi2_subchannels(pars); };

        // Chi2 contribution from just saved sub-channel data
        double chi2_subchannels(const double *par);

        // --------------------------------------------------------------------
        // Parameter handling

        int _N_pars = 0;   // total # of parameters
        int _N_amp  = 0;   // # parameters in production line-shape
        int _N_V    = 0;   // # parameters in decay amplitude

        // Container struct carrying all relevant info for each parameter involved in the fit
        struct parameter
        {
            parameter(int i)
            : _i(i), _label("par["+std::to_string(i)+"]")
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

        // Input a variable label and find its corresponding index
        int find_parameter(std::string name)
        {
            for (int i = 0; i < _pars.size(); i++)
            {
                if (_pars[i]._label != name) return i;
            }

            warning("amplitude_fitter::find_parameter()", "Parameter named " + name + " not found!");
            return -1;
        };

        // convert double[] to vector<double>'s
        std::vector<double> convert(const double * par);

        // This converts but also splits the single vector into two depending on the parameters which go into _V and which go to _amp
        // And passes the appropriate length vectors to amplitudes with amplitude::set_parameters
        void allocate_parameters(const double *par);

        // --------------------------------------------------------------------
        // Data handling

        // Running total number of data points saved
        int _N_data = 0;

        // Simple container struct to hold all the relevent info for each user-added data set to fit against
        struct data_set
        {
            data_set(subchannel abc, double sqrts, std::vector<double> sqrtsigmas, std::vector<double> data, std::array<std::vector<double>,2> errors, std::string id)
            : _subchannel(abc), _sqrts(sqrts), _sqrtsigmas(sqrtsigmas), _data(data), _errors(errors), _id(id),
            _N(sqrtsigmas.size())
            {};

            std::string _id;
            subchannel _subchannel;
            double _sqrts;
            std::vector<double> _sqrtsigmas, _data;
            std::array<std::vector<double>,2> _errors;
            int _N;
            double _norm;
        };
        // Vector containing all our different data sets
        std::vector<data_set> _subchannel_data;

        // --------------------------------------------------------------------
        // Status message handling

        // Methods for printing out messages
        inline void divider()
        {
            std::cout << "--------------------------------------------------------------" << std::endl;
        };
        inline void new_line()
        {
            std::cout << std::endl;
        };

        // Simply return string for default particle labels
        inline std::string subchannel_label(subchannel x)
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
};

#endif