// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude_fitter.hpp"

// ---------------------------------------------------------------------------
// Primary function, sets up the fitter with previously set parameters and data sets
// Starts fit with an user-supplied vector of starting values for all parameters
void hadMolee::amplitude_fitter::do_fit(std::vector<double> starting_guess)
{  
    if (starting_guess.size() != _pars.size())
    {
        warning("amplitude_fitter::do_fit()", "Size of initial guess vector doesnt match number of parameters!");
        return;
    };

    // If size of vector matches, feed all parameter info to minuit
    set_up(starting_guess);

    // Print out info on variables and data to command line
    new_line(); data_info();
    new_line(); divider(); 
    variable_info(starting_guess, 0);
    new_line(); divider(); new_line();

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Beginning fit..." << std::flush; 

    if (_printLevel != 0) new_line();   
    _minuit->Minimize();
    if (_printLevel != 0) new_line();   

    std::cout << "Done! \n";
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
    std::cout << std::left << "Finished in " << duration.count() << " sec" << std::endl;

    print_results();

    return;
};

// Minimization function for the subchannel data
double hadMolee::amplitude_fitter::chi2_subchannels(const double *par)
{   
    allocate_parameters(par);

    // Iterate over each saved dataset and add to global chi2 
    double chi2 = 0.;
    for (data_set data : _subchannel_data)
    {
        // Fixed e+e- invariant mass
        double s     = pow(data._sqrts, 2.);

        // We assume the subchannel data set is in raw event counts
        // Thus we normalize our curves to have the same area        
        _amplitude->normalize(data._norm, s);
        
        // Loop over data points calculating chi2
        double chi2_i = 0.;
        for (int i = 0; i < data._N; i++)
        {
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
void hadMolee::amplitude_fitter::add_subchannel_data(subchannel abc, double sqs, std::vector<double> sqsig, std::vector<double> data, std::array<std::vector<double>,2> errors, std::string id)
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

    // Calculate normalization
    double sum = 0;
    for (double data : data){ sum += data; };
    new_data._norm = sum;

    _subchannel_data.push_back(new_data);

    // Add number of points to the running totals
    _N_data += n;
};

// Alternatively you can jsut specify the subchannel, fixed center-of-mass energy, and point to a data file
void hadMolee::amplitude_fitter::add_subchannel_data(subchannel abc, double sqs, std::string filename, std::string id)
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
};

//-----------------------------------------------------------------------
// Parameter management methods 

void hadMolee::amplitude_fitter::set_parameter_labels(std::vector<std::string> labels)
{
    if (labels.size() != _pars.size())
    {
        std::cout << "amplitude_fitter::set_parameter_labels() : Input vector is of incorrect size! Expected " << _pars.size() << " but recieved " << labels.size() << "!" << std::endl;
        std::cout << "Continuing without adding labels..." << std::endl;
        return;
    };

    for (int i = 0; i < _pars.size(); i++)
    {
        _pars[i]._label = labels[i];
    }
};

void hadMolee::amplitude_fitter::set_parameter_limits(int i, std::array<double,2> ranges, double step)
{
    _pars[i]._custom_limits = true;
    _pars[i]._lower_limit   = ranges[0];
    _pars[i]._upper_limit   = ranges[1];
    _pars[i]._step_size     = step;
};

// convert double[] to vector<double>'s
// Also splits the single vector into two depending on the parameters which go into _V and which go to _amp
// This splitting is necessary in case multiple _amps all share the same _V
std::vector<double> hadMolee::amplitude_fitter::convert(const double * par)
{
    std::vector<double> all_pars;
    for (int n = 0; n < _N_pars; n++)
    {
        all_pars.push_back(par[n]);
    };

    return all_pars;
};

// Takes in the double* from minuit, and feeds the vectors of appropriate size to each subamplitude
void hadMolee::amplitude_fitter::allocate_parameters(const double *par)
{
    // Split the parameter vector 
    std::vector<double> all_pars = convert(par);

    // now we split into two sets
    // The first are the parameters that belong to our production lineshape model
    auto start = all_pars.begin();
    auto end   = all_pars.begin() + _N_V;

    std::vector<double> V_pars(start, end);

    // The second are for the decay amplitude
    start = all_pars.begin() + _N_V;
    end   = all_pars.end();

    std::vector<double> amp_pars(start, end);

    // And allocate them appropriately
    _V->set_parameters(V_pars); 
    _amplitude->set_parameters(amp_pars);
};

//-----------------------------------------------------------------------
// Methods to print out relevant info to command line

// Print out a summary of saved data
void hadMolee::amplitude_fitter::data_info()
{
    divider();
    std::cout << std::left << std::setw(30) << "Using e+e- lineshape model:" << std::setw(20) << _V->get_id() << std::endl;
    new_line();

    std::cout << std::left << "Fitting decay amplitude for " << _amplitude->_kinematics->get_id() << " final-state (\"" << _amplitude->get_id() << "\") to " << _N_data << " data points: \n";
    new_line();
    std::cout << std::left << std::setw(25) << "DATA SET"            << std::setw(20) << "CHANNEL"         << std::setw(10) << "POINTS"  << std::endl;
    std::cout << std::left << std::setw(25) << "----------------"    << std::setw(20) << "--------------"  << std::setw(10) << "-------" << std::endl;

    for (int k = 0; k < _subchannel_data.size(); k++)
    {
        std::cout << std::left << std::setw(25) << _subchannel_data[k]._id  
                               << std::setw(20) << _amplitude->_kinematics->subchannel_label(_subchannel_data[k]._subchannel)  
                               << std::setw(10) << _subchannel_data[k]._data.size()  << std::endl;  
    };
};

void hadMolee::amplitude_fitter::set_up(std::vector<double> starting_guess)
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

    fcn = ROOT::Math::Functor(this, &amplitude_fitter::chi2, _pars.size());
    _minuit->SetFunction(fcn);
};


// Print out a little table of the current status of parameters
void hadMolee::amplitude_fitter::variable_info(std::vector<double> starting_guess, bool opt)
{  
    std::cout << std::setprecision(10);

    std::string column_3;
    (opt) ? (column_3 = "FIT VALUE") : (column_3 = "START VALUE");

    if (!opt)  std::cout << std::left << "Fitting " + std::to_string(_minuit->NFree()) << " (out of " << std::to_string(_minuit->NDim()) << ") parameters:\n" << std::endl; 
    std::cout << std::left << std::setw(10) << "N"       << std::setw(20) << "PARAMETER"  << std::setw(10) << column_3       << std::endl;
    std::cout << std::left << std::setw(10) << "-----"   << std::setw(20) << "----------" << std::setw(10) << "------------" << std::endl;

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
        std::cout << std::left << std::setw(10) << i << std::setw(20) << _pars[i]._label << std::setw(20) << starting_guess[i] << std::setw(10) << extra << std::endl;
    };
};

void hadMolee::amplitude_fitter::print_results()
{
    // Fit results
    int dof = _N_pars - _minuit->NFree();

    new_line();
    double chi2    = _minuit->MinValue();
    double chi2dof = chi2 / double(dof) ;
    std::vector<double> best_params = convert(_minuit->X());

    divider();
    std::cout << std::left << std::setw(10) << "chi2 = "      << std::setw(15) << chi2;
    std::cout << std::left << std::setw(10) << "chi2/dof = "  << std::setw(15) << chi2dof << std::endl;
    new_line();
    variable_info(best_params, 1);
    divider();
    new_line();
    
    // At the end update the amplitude parameters to include the fit results
    allocate_parameters(_minuit->X());
}