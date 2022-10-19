// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "amplitude_fitter.hpp"

void amplitude_fitter::do_fit()
{  
    data_info();
};

//-----------------------------------------------------------------------
// Methods to print out relevant info to command line

// Print out a summary of saved data
void amplitude_fitter::data_info()
{
    cout << left << "Fitting amplitude (\"" << _amplitude->get_id() << "\") to " << _N << " data points: \n";
    new_line();
    cout << left << setw(30) << "DATA SET"            << setw(20) << "SUBCHANNEL"      << setw(10) << "POINTS"  << endl;
    cout << left << setw(30) << "----------------"    << setw(20) << "--------------"  << setw(10) << "-------" << endl;

    for (int k = 0; k < _subchannel_data.size(); k++)
    {
        cout << left << setw(30) << _subchannel_data[k]._id  << setw(20)  << subchannel_label(_subchannel_data[k]._subchannel)   << setw(10) <<  _subchannel_data[k]._data.size()  << endl;  
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
    int dof = _N - _minuit->NFree();

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
    _amplitude->set_parameters(best_params);
}