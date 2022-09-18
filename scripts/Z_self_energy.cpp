#include "hadronic_molecule.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void Z_self_energy()
{   
    // Dressed propagator
    DsD_molecule Z_meson;
    
    // ------------------------------------------------------
    double W = M_Y;
    double s = W*W;

    int N = 100;
    bool PRINT = true;

    double xmin = 3.8;
    double xmax = 4.1;

    double ymin = -1.;
    double ymax = 2.;

    string xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";
    string ylabel   = "#it{i}#Sigma_{#it{D}*#it{D}} #times 100";
    string filename = "self_energy.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    auto reSigma = [&] (double E)
    {
        return 100*real(Z_meson.self_energy(E));
    };
    auto imSigma = [&] (double E)
    {
        return 100*imag(Z_meson.self_energy(E));
    };
    plotter->AddEntry(N, reSigma, {xmin, xmax}, "Real", PRINT);
    plotter->AddEntry(N, imSigma, {xmin, xmax}, "Imaginary", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.25,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;

};  