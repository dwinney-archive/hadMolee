
#include "reaction_kinematics.hpp"
#include "DsDpi.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void tree_widths()
{   
    
    double W = 4.260;
    double s = W*W;

    reaction_kinematics kDsDpi({M_DSTAR, M_D, M_PION});
    kDsDpi.set_particle_labels({"#it{D}*", "#it{D}", "#pi"});

    // Empty Y propagator
    hadronic_molecule Y_meson;

    DsDpi_tree tree(&kDsDpi, &Y_meson);

    int N = 100;
    bool PRINT = false;

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D(0);

    string xlabel, ylabel, filename;
    double xmin, xmax, ymin, ymax;

    // ---------------------------------------------------------------------------
    // Print everything out

    auto dgammaab = [&](double w)
    {
        double sigma = w * w;
        return tree.dGamma_ab(s, sigma);
    };

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_b());
    xmax = (W - kDsDpi.mass_c());

    ymin = 0.;
    ymax = 4.;  

    xlabel = "#sqrt{#sigma_{"+ kDsDpi.particle_a() + kDsDpi.particle_b() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}  [a.u.]";
    filename = "tree_width_ab.pdf";

    tree.set_debug(0); tree.normalize(1., s);
    plotter->AddEntry(N, dgammaab, {xmin, xmax}, "S and D waves", PRINT);
    tree.set_debug(1); tree.normalize(1., s);
    plotter->AddEntry(N, dgammaab, {xmin, xmax}, "D wave only", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.6,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    plotter->ClearData();

    auto dgammabc = [&](double w)
    {
        double sigma = w * w;
        return tree.dGamma_bc(s, sigma);
    };

    xmin = (kDsDpi.mass_b() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_a());

    ymax = 4.;

    xlabel = "#sqrt{#sigma_{"+ kDsDpi.particle_b() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_b() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "tree_width_bc.pdf";

    tree.set_debug(0); tree.normalize(1., s);
    plotter->AddEntry(N, dgammabc, {xmin, xmax}, "S and D waves", PRINT);
    tree.set_debug(1); tree.normalize(1., s);
    plotter->AddEntry(N, dgammabc, {xmin, xmax}, "D wave only", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.3,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    plotter->ClearData();

    auto dgammaac = [&](double w)
    {
        double sigma = w * w;
        return tree.dGamma_ac(s, sigma);
    };

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_b());

    ymax = 7.;

    xlabel = "#sqrt{#sigma_{"+ kDsDpi.particle_a() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "tree_width_ac.pdf";

    tree.set_debug(0); tree.normalize(1., s);
    plotter->AddEntry(N, dgammaac, {xmin, xmax}, "S and D waves", PRINT);
    tree.set_debug(1); tree.normalize(1., s);
    plotter->AddEntry(N, dgammaac, {xmin, xmax}, "D wave only", PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.3,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    delete plotter;

};  