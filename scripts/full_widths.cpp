
#include "reaction_kinematics.hpp"
#include "amplitude_sum.hpp"
#include "DsDpi.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void full_widths()
{   
    // ---------------------------------------------------------------------------
    // Amplitude definition and set up

    // Fix the e+e- mass to the nominal Y mass
    double W = M_Y4260;
    double s = W*W;

    // Define the global kinematics instance for the D*Dpi final state
    reaction_kinematics kDsDpi({M_DSTAR, M_D, M_PION});
    kDsDpi.set_particle_labels({"#it{D}*", "#it{D}", "#pi"});

    // Empty Y propagator since we only care about the widths at fixed W
    D1D_molecule Y;

    // S-wave we have a single term

    DsDpi_swave swave(&kDsDpi, &Y);
    swave.set_parameters({-12.67, -15.23});

    // D-wave recieves contributions from two diagrams: 

    // the triangle
    DsDpi_triangle triangle(&kDsDpi, &Y, "Triangle");
    triangle.set_debug(1); // Only D-wave couplings

    // and the D1 tree transition
    DsDpi_tree tree(&kDsDpi, &Y, "Tree");
    tree.set_debug(1); // Only D-wave couplings

    // These can be summed up like this:
    amplitude_sum full(&kDsDpi, &Y, "D-wave");
    full.add({&tree, &triangle, &swave});
    // full.normalize(1., s);

    // ---------------------------------------------------------------------------
    // Plotter object and options

    jpacGraph1D * plotter = new jpacGraph1D(0);

    int N = 100;
    bool PRINT = true;
    string xlabel, ylabel, filename;
    double xmin, xmax, ymin, ymax;

    auto dgammaab = [&](double w)
    {
        double sigma = w * w;
        return full.dGamma_ab(s, sigma);
    };

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_b());
    xmax = (W - kDsDpi.mass_c());

    ymin = 0.;
    ymax = 7.;  

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}  [a.u.]";
    filename = "full_width_ab.pdf";

    plotter->AddEntry(N, dgammaab, {xmin, xmax}, "D-wave", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(false);
    plotter->Plot(filename);

    plotter->ClearData();

    auto dgammabc = [&](double w)
    {
        double sigma = w * w;
        return full.dGamma_bc(s, sigma);
    };

    xmin = (kDsDpi.mass_b() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_a());

    ymax = 7.;

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_b() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_b() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "full_width_bc.pdf";

    plotter->AddEntry(N, dgammabc, {xmin, xmax}, "D-wave", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(0.3,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    plotter->ClearData();

    auto dgammaac = [&](double w)
    {
        double sigma = w * w;
        return full.dGamma_ac(s, sigma);
    };

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_b());

    ymax = 7.;

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_a() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "full_width_ac.pdf";

    plotter->AddEntry(N, dgammaac, {xmin, xmax}, "D-wave", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(0.3,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    delete plotter;
};  