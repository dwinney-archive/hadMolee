#include "hadronic_molecule.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void Y_self_energy()
{   
    // Dressed propagator
    D1D_molecule Y;
    
    // ------------------------------------------------------

    int N = 100;
    bool PRINT = true;

    double xmin = 4.1;
    double xmax = 4.4;

    double ymin = -1.5;
    double ymax = 2.;

    string xlabel   = "#sqrt{s}  [GeV]";
    string ylabel   = "#it{i}#Delta_{Y}";
    string filename = "self_energy.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    Y.set_pole_mass(4.26);

    auto f = [&] (double x)
    {
        return Y.propagator(x);
    };

    // Split the above function above into re and im
    auto reSigma = [&] (double E)
    {
        return real(f(E));
    };
    auto imSigma = [&] (double E)
    {
        return imag(f(E));
    };
    plotter->AddEntry(N, reSigma, {xmin, xmax}, "Real", PRINT);
    plotter->AddEntry(N, imSigma, {xmin, xmax}, "Imaginary", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.25,0.2);
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;

};  