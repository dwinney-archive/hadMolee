// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "fitter.hpp"

namespace hadMolee
{
    // ---------------------------------------------------------------------------
    // Primary function, sets up the fitter with previously set parameters and data sets
    // Starts fit with an user-supplied vector of starting values for all parameters
    void fitter::do_fit(std::vector<double> starting_guess)
    {  
        if (_N_data == 0)
        {
            warning("fitter::do_fit()", "No data-points saved! Returning without fit...");
            return;
        };

        if (starting_guess.size() != _pars.size())
        {
            warning("fitter::do_fit()", "Size of initial guess vector doesnt match number of parameters!");
            return;
        };

        // If size of vector matches, feed all parameter info to minuit
        set_up(starting_guess);

        // Print out info on variables and data to command line
        empty_line(); data_info();
        empty_line(); divider(); 
        variable_info(starting_guess, 0);
        empty_line(); divider(); empty_line();
    
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Beginning fit..." << std::flush; 

        if (_printLevel != 0) empty_line();   
        _minuit->Minimize();
        if (_printLevel != 0) empty_line();   

        std::cout << "Done! \n";
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
        std::cout << std::left << "Finished in " << duration.count() << " sec" << std::endl;

        print_results();

        return;
    };

    // Minimization function for the subchannel data
    double fitter::chi2_subchannels(const double *par)
    {   
        allocate_parameters(par);

        // Iterate over each saved dataset and add to global chi2 
        double chi2 = 0.;
        for (data_set data : _subchannel_data)
        {
            // Fixed e+e- invariant mass
            double s     = pow(data._sqrts, 2.);

            // Loop over data points calculating chi2
            double chi2_i = 0.;
            for (int i = 0; i < data._N; i++)
            {
                // Apply the subchannel normalization to the amplitude
                _amplitude->normalize(data._normalization);

                // Fixed sub-channel invariant mass
                double sigma = pow(data._sqrtsigmas[i], 2.);

                // Caculate intensity from saved amplitude and from data
                double I_th   = _amplitude->dGamma(data._subchannel, s, sigma);
                double I_ex   = data._data[i];
                double error  = (data._errors[0][i] + data._errors[1][i]);

                // Add to local chi2
                chi2_i += pow( (I_th - I_ex) / error, 2.);
            };

            // Add the local chi2 to the running chi2
            chi2 += chi2_i;
        };

        return chi2;
    };

    //-----------------------------------------------------------------------
    // Data management methods 

    // Add vectors corresponding to partial-width data in a specific subchannel
    // Arguments require:
    // - specify subchannel according to labels specified in amp->_reaction_kinematics
    // - fixed sqrt(s) of the photon 
    // - vector of sqrt(sigma) where sigma is invairant mass of subsystem
    // - vector data points
    // - vector of errors
    // String id optional parameter for feeding fit results to a plotter object
    void fitter::add_subchannel_data(subchannel abc, double sqs, std::vector<double> sqsig, std::vector<double> data, std::array<std::vector<double>,2> errors, std::string id)
    {
        // Number of data points for this set
        int n = sqsig.size();

        if (data.size() != n || errors[0].size() != n || errors[1].size() != n) 
        {
            warning("fitter::add_subchannel_data", "Vectors received not the correct size! Skipping data set " + id + "...");
            return;
        };

        if ( id == "" ) id = subchannel_label(abc) + "_data[" + std::to_string(_subchannel_data.size()) + "]";
        data_set new_data(_subchannel_data.size(), abc, sqs, sqsig, data, errors, id);

        _subchannel_data.push_back(new_data);

        // Add number of points to the running totals
        _N_data += n;
    };

    // Alternatively you can jsut specify the subchannel, fixed center-of-mass energy, and point to a data file
    void fitter::add_subchannel_data(subchannel abc, double sqs, std::string filename, std::string id)
    {
        // Import data
        // The resulting format will be in 4 vectors
        auto file = import_data(filename);

        // The first two columns are already in the format we need
        std::vector<double> sqsig  = file[0];
        std::vector<double> data   = file[1];
        std::array<std::vector<double>,2> errors = {file[2], file[3]};
        
        // Save to fitter
        add_subchannel_data(abc, sqs, sqsig, data, errors, id);    

        // Increment count of normalization parameters needed
        _N_norms++;
    };

    //-----------------------------------------------------------------------
    // Parameter management methods 

    void fitter::set_parameter_labels(std::vector<std::string> labels)
    {
        if (labels.size() != _pars.size())
        {
            std::cout << "fitter::set_parameter_labels() : Input vector is of incorrect size! Expected " << _pars.size() << " but recieved " << labels.size() << "!" << std::endl;
            std::cout << "Continuing without adding labels..." << std::endl;
            return;
        };

        for (int i = 0; i < _pars.size(); i++)
        {
            _pars[i]._label = labels[i];
        }
    };

    void fitter::set_parameter_limits(int i, std::array<double,2> ranges, double step)
    {
        _pars[i]._custom_limits = true;
        _pars[i]._lower_limit   = ranges[0];
        _pars[i]._upper_limit   = ranges[1];
        _pars[i]._step_size     = step;
    };

    // convert double[] to vector<double>'s
    // Also splits the single vector into two depending on the parameters which go into _V and which go to _amp
    // This splitting is necessary in case multiple _amps all share the same _V
    std::vector<double> fitter::convert(const double * par)
    {
        std::vector<double> all_pars;
        for (int n = 0; n < _N_pars + _N_norms; n++)
        {
            all_pars.push_back(par[n]);
        };

        return all_pars;
    };

    // Takes in the double* from minuit, and feeds the vectors of appropriate size to each subamplitude
    void fitter::allocate_parameters(const double *par, bool best_fit)
    {
        // Split the parameter vector 
        std::vector<double> all_pars = convert(par);

        // now we split into three sets
        // The first are the parameters that belong to our production lineshape model
        auto start = all_pars.begin();
        auto end   = all_pars.begin() + _N_V;

        std::vector<double> V_pars(start, end);

        // The second are for the decay amplitude
        start = all_pars.begin() + _N_V;
        end   = start + _N_amp;

        std::vector<double> amp_pars(start, end);

        // Lastly are the data-set normalizations
        start = all_pars.begin() + _N_V + _N_amp;
        end   = start + _N_norms;
        std::vector<double> norms(start, end);

        // And allocate them appropriately
        _V->set_parameters(V_pars); 
        _amplitude->set_parameters(amp_pars);
        for (int i = 0; i < _N_norms; i++)
        {
            _subchannel_data[i]._normalization = norms[i];
        };

        // The bool best_fit saves the vector to member data so it can be output
        if (best_fit) 
        {
            start = all_pars.begin();
            end   = all_pars.begin() + _N_V + _N_amp;
            _best_fit = std::vector<double>(start, end);

            _normalizations = norms;

            // In case some normalization of the amplitude is lingering we reset it 
            _amplitude->normalize(1.);
        }
    };

    //-----------------------------------------------------------------------
    // Methods to print out relevant info to command line

    // Print out a summary of saved data
    void fitter::data_info()
    {
        divider();
        std::cout << std::left << std::setw(28) << "Using e+e- lineshape model:" << std::setw(_V->get_id().length() + 2) << _V->get_id() << " (" + std::to_string(_V->N_parameters()) + " free parameters)" << std::endl;
        empty_line();
        std::cout << std::left << "Using decays to 1 final-state with total of " << _N_data << " data points: \n";
        empty_line();
        dashed_divider();
        centered("FINAL-STATE: " +  _amplitude->_kinematics->get_id());
        empty_line();
        std::cout << std::left << std::setw(13) <<  "AMPLITUDE: " << std::setw(_amplitude->get_id().length() + 2) << _amplitude->get_id() << " (" + std::to_string(_amplitude->N_parameters()) + " free parameters)" << std::endl;
        empty_line();
        
        std::cout << std::left << std::setw(23) << "DATA SET"               << std::setw(11) << "SQRT(s)"  << std::setw(13) << "CHANNEL"     << std::setw(10) << "POINTS"  << std::endl;
        std::cout << std::left << std::setw(23) << "-------------------"   << std::setw(11) << "--------" << std::setw(13) << "----------"  << std::setw(10) << "-------" << std::endl;

        for (data_set dataset : _subchannel_data)
        {
            std::cout << std::left  << std::setw(23) << dataset._id  
                                    << std::setw(11) << dataset._sqrts
                                    << std::setw(13) << _amplitude->_kinematics->subchannel_label(dataset._subchannel)  
                                    << std::setw(10) << dataset._data.size()  << std::endl;  
        };
    };

    void fitter::set_up(std::vector<double> starting_guess)
    {
        _minuit->Clear();
        _minuit->SetTolerance(_tolerance);
        _minuit->SetPrintLevel(_printLevel);
        _minuit->SetMaxFunctionCalls(_maxCalls);

        for (int a = 0; a < _pars.size(); a++)
        {   
            _minuit->SetVariable(a, _pars[a]._label, starting_guess[a], _pars[a]._step_size);

            if (_pars[a]._custom_limits)
            {
                _minuit->SetVariableLimits(a, _pars[a]._lower_limit, _pars[a]._upper_limit);
            }
            if (_pars[a]._fixed)
            {
                _minuit->FixVariable(a);
            }
        };
        
        // For each subchannel dataset we need to include one normalization as a fitting parameter
        for (int b = 0; b < _N_norms; b++)
        {
            // Give each N a default name and start at 1
            std::string norm_label = "N["+ std::to_string(b) + "]";
            _minuit->SetVariable(_pars.size() + b, norm_label, 1., 0.1);

        };
        
        fcn = ROOT::Math::Functor(this, &fitter::chi2, _pars.size() + _subchannel_data.size());
        _minuit->SetFunction(fcn);
    };


    // Print out a little table of the current status of parameters
    void fitter::variable_info(std::vector<double> starting_guess, bool opt)
    {  
        std::cout << std::setprecision(10);

        std::string column_3;
        (opt) ? (column_3 = "FIT VALUE") : (column_3 = "START VALUE");

        if (!opt)  std::cout << std::left << "Fitting " + std::to_string(_minuit->NFree() - _N_norms) << " (out of " << std::to_string(_minuit->NDim() - _N_norms) << ") model parameters and " << _subchannel_data.size() << " normalizations:" << std::endl;
        empty_line();

        std::cout << std::left << std::setw(7) << "N"       << std::setw(20) << "PARAMETER"  << std::setw(10) << column_3       << std::endl;
        std::cout << std::left << std::setw(7) << "----"   << std::setw(20) << "----------" << std::setw(10) << "------------" << std::endl;

        for (int i = 0; i < _pars.size(); i++)
        {
            std::string extra = "";
            if (_pars[i]._custom_limits && !opt)
            {   
                std::stringstream ss;
                ss << std::setprecision(5) << "[" << _pars[i]._lower_limit << ", " << _pars[i]._upper_limit << "]";
                extra = ss.str();
            };
            if (_pars[i]._fixed && !opt) extra = "FIXED";
            std::cout << std::left << std::setw(7) << i << std::setw(20) << _pars[i]._label << std::setw(20) << starting_guess[i] << std::setw(10) << extra << std::endl;
        };

        if (opt)
        {
            empty_line();
            std::cout << std::left << std::setw(23) << "DATA SET"             << std::setw(11) << "SQRT(s)"  << std::setw(13) << "CHANNEL"     << std::setw(15) << "NORMALIZATION"  << std::endl;
            std::cout << std::left << std::setw(23) << "-------------------"   << std::setw(11) << "--------" << std::setw(13) << "----------"  << std::setw(15) << "--------------" << std::endl;

            for (data_set dataset : _subchannel_data)
            {
                std::cout << std::left  << std::setw(23) << dataset._id  
                                        << std::setw(11) << dataset._sqrts
                                        << std::setw(13) << _amplitude->_kinematics->subchannel_label(dataset._subchannel)  
                                        << std::setw(15) << dataset._normalization << std::endl;  
            };
        }
    };

    void fitter::print_results()
    {
        // Fit results
        int dof = _N_data - _minuit->NFree();

        empty_line();
        double chi2    = _minuit->MinValue();
        double chi2dof = chi2 / double(dof) ;
        std::vector<double> best_params = convert(_minuit->X());

        divider();
        std::cout << std::left << std::setw(10) << "chi2 = "      << std::setw(15) << chi2;
        std::cout << std::left << std::setw(10) << "chi2/dof = "  << std::setw(15) << chi2dof << std::endl;
        empty_line();
        variable_info(best_params, 1);
        empty_line();
        divider();
        empty_line();
        
        // At the end update the amplitude parameters to include the fit results
        allocate_parameters(_minuit->X(), true);
    }
};