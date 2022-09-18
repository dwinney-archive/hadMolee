#include "hadronic_molecule.hpp"
#include "propagator.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void Z_propagator()
{   
    // Dressed propagator
    DsD_molecule Z_meson;

    // Nonrelativistic BW for comparison
    nonrelativistic_BW Z_BW(M_Z, 28.4E-3);
    
    // ------------------------------------------------------
    double W = M_Y;
    double s = W*W;

    int N = 200;
    bool PRINT = true;

    double xmin = 3.8;
    double xmax = 4.1;

    double ymin = -20.;
    double ymax =  40.;

    string xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";
    string ylabel   = "#it{i}#it{G}_{#it{Z}}";
    string filename = "Z_prop.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    auto reiG = [&] (double E)
    {
        return real(XI*Z_meson.propagator(E));
    };
    auto imiG = [&] (double E)
    {
        return imag(XI*Z_meson.propagator(E));
    };

    auto reiBW = [&] (double E)
    {
        return real(XI*Z_BW.eval(E));
    };
    auto imiBW = [&] (double E)
    {
        return imag(XI*Z_BW.eval(E));
    };

    plotter->AddEntry(N, reiG, {xmin, xmax}, "Real", PRINT);
    plotter->AddDashedEntry(N, reiBW, {xmin, xmax},  PRINT);

    plotter->AddEntry(N, imiG, {xmin, xmax}, "Imaginary", PRINT);
    plotter->AddDashedEntry(N, imiBW, {xmin, xmax},  PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.6,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;
};  