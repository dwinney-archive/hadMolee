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
void amplitude_fitter::do_fit(vector<double> starting_guess)
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

    auto start = chrono::high_resolution_clock::now();
    cout << "Beginning fit..." << flush; 

    if (_printLevel != 0) new_line();   
    _minuit->Minimize();
    if (_printLevel != 0) new_line();   

     cout << "Done! \n";
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast< chrono::seconds>(stop - start);
    cout << left << "Finished in " << duration.count() << " sec" << endl;

    print_results();

    return;
};

// Minimization function for the subchannel data
double amplitude_fitter::chi2_subchannels(const double *par)
{   
    allocate_parameters(par);

    // Iterate over each saved dataset and add to global chi2 
    double chi2 = 0.;
    for (data_set data : _subchannel_data)
    {
        // We assume the subchannel data set is in raw event counts
        // Thus we normalize our curves to have the same area
        double total_events = 0;
        for (auto bin : data._data){ total_events += bin; }
        
        // Fixed e+e- invariant mass
        double s     = pow(data._sqrts, 2.);

        _amplitude->normalize(total_events, s);
        
        // Loop over data points calculating chi2
        double chi2_i = 0.;
        for (int i = 0; i < data._N; i++)
        {
            // Fixed sub-channel invariant mass
            double sigma = pow(data._sqrtsigmas[i], 2.);

            // Caculate intensity from saved amplitude and from data
            double I_th   = _amplitude->dGamma(data._subchannel, s, sigma);
            double I_ex   = data._data[i];
            double error  = data._errors[i];

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
void amplitude_fitter::add_subchannel_data(subchannel abc, double sqs, vector<double> sqsig, vector<double> data, vector<double> errors, string id)
{
    // Number of data points for this set
    int n = sqsig.size();

    if (data.size() != n || errors.size() != n) 
    {
        warning("amplitude_fitter::add_subchannel_data", "Vectors received not the correct size! Skipping data set " + id + "...");
        return;
    };

    if ( id == "" ) id = subchannel_label(abc) + "_data[" + to_string(_subchannel_data.size()) + "]";
    data_set new_data(abc, sqs, sqsig, data, errors, id);

    _subchannel_data.push_back(new_data);

    // Add number of points to the running totals
    _N_data += n;
};

// Alternatively you can jsut specify the subchannel, fixed center-of-mass energy, and point to a data file
void amplitude_fitter::add_subchannel_data(subchannel abc, double sqs, string filename, string id)
{
    // Import data
    // The resulting format will be in 4 vectors
    auto file = import_data(filename);

    // The first two columns are already in the format we need
    std::vector<double> sqsig = file[0];
    std::vector<double> data  = file[1];

    // The errors are given as upper and lower values so we neeed to convert them to a single |error| = |upper - lower|
    std::vector<double>  errors;
    for (int i = 0; i < file[2].size(); i++)
    {
        errors.push_back( file[3][i] - file[2][i] );
    };

    // Save to fitter
    add_subchannel_data(abc, sqs, sqsig, data, errors, id);    
};

//-----------------------------------------------------------------------
// Parameter management methods 

void amplitude_fitter::set_parameter_labels(vector<string> labels)
{
    if (labels.size() != _pars.size())
    {
        cout << "amplitude_fitter::set_parameter_labels() : Input vector is of incorrect size! Expected " << _pars.size() << " but recieved " << labels.size() << "!" << endl;
        cout << "Continuing without adding labels..." << endl;
        return;
    };

    for (int i = 0; i < _pars.size(); i++)
    {
        _pars[i]._label = labels[i];
    }
};

void amplitude_fitter::set_parameter_limits(int i, array<double,2> ranges, double step)
{
    _pars[i]._custom_limits = true;
    _pars[i]._lower_limit   = ranges[0];
    _pars[i]._upper_limit   = ranges[1];
    _pars[i]._step_size     = step;
};

// convert double[] to vector<double>'s
// Also splits the single vector into two depending on the parameters which go into _V and which go to _amp
// This splitting is necessary in case multiple _amps all share the same _V
vector<double> amplitude_fitter::convert(const double * par)
{
    vector<double> all_pars;
    for (int n = 0; n < _N_pars; n++)
    {
        all_pars.push_back(par[n]);
    };

    return all_pars;
};

// Takes in the double* from minuit, and feeds the vectors of appropriate size to each subamplitude
void amplitude_fitter::allocate_parameters(const double *par)
{
    // Split the parameter vector 
    std::vector<double> all_pars = convert(par);

    // now we split into two sets
    // The first are the parameters that belong to our production lineshape model
    auto start = all_pars.begin();
    auto end   = all_pars.begin() + _N_V;

    vector<double> V_pars(start, end);

    // The second are for the decay amplitude
    start = all_pars.begin() + _N_V;
    end   = all_pars.end();

    vector<double> amp_pars(start, end);

    // And allocate them appropriately
    _V->set_parameters(V_pars); 
    _amplitude->set_parameters(amp_pars);
};

//-----------------------------------------------------------------------
// Methods to print out relevant info to command line

// Print out a summary of saved data
void amplitude_fitter::data_info()
{
    divider();
    cout << left << setw(30) << "Using e+e- lineshape model:" << setw(20) << _V->get_id() << endl;
    new_line();

    cout << left << "Fitting decay amplitude for " << _amplitude->_kinematics->get_id() << " final-state (\"" << _amplitude->get_id() << "\") to " << _N_data << " data points: \n";
    new_line();
    cout << left << setw(25) << "DATA SET"            << setw(20) << "CHANNEL"     << setw(10) << "POINTS"  << endl;
    cout << left << setw(25) << "----------------"    << setw(20) << "--------------"  << setw(10) << "-------" << endl;

    for (int k = 0; k < _subchannel_data.size(); k++)
    {
        cout << left << setw(25) << _subchannel_data[k]._id  
                     << setw(20) << _amplitude->_kinematics->subchannel_label(_subchannel_data[k]._subchannel)  
                     << setw(10) << _subchannel_data[k]._data.size()  << endl;  
    };
};

void amplitude_fitter::set_up(vector<double> starting_guess)
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
void amplitude_fitter::variable_info(vector<double> starting_guess, bool opt)
{  
    cout << setprecision(10);

    string column_3;
    (opt) ? (column_3 = "FIT VALUE") : (column_3 = "START VALUE");

    if (!opt)  cout << left << "Fitting " + to_string(_minuit->NFree()) << " (out of " << to_string(_minuit->NDim()) << ") parameters:\n" << endl; 
    cout << left << setw(10) << "N" << setw(20) << "PARAMETER" << setw(10) << column_3 << endl;
    cout << left << setw(10) << "-----" << setw(20) << "----------" << setw(10) << "------------"<< endl;

    for (int i = 0; i < _pars.size(); i++)
    {
        string extra = "";
        if (_pars[i]._custom_limits && !opt)
        {   
            stringstream ss;
            ss << setprecision(5) << "[" << _pars[i]._lower_limit << ", " << _pars[i]._upper_limit << "]";
            extra = ss.str();
        };
        if (_pars[i]._fixed && !opt) extra = "FIXED";
        cout << left << setw(10) << i << setw(20) << _pars[i]._label << setw(20) << starting_guess[i] << setw(10) << extra << endl;
    };
};

void amplitude_fitter::print_results()
{
    // Fit results
    int dof = _N_pars - _minuit->NFree();

    new_line();
    double chi2    = _minuit->MinValue();
    double chi2dof = chi2 / double(dof) ;
    vector<double> best_params = convert(_minuit->X());

    divider();
    cout << left << setw(10) << "chi2 = "      << setw(15) << chi2;
    cout << left << setw(10) << "chi2/dof = "  << setw(15) << chi2dof << endl;
    new_line();
    variable_info(best_params, 1);
    divider();
    new_line();
    
    // At the end update the amplitude parameters to include the fit results
    allocate_parameters(_minuit->X());
}