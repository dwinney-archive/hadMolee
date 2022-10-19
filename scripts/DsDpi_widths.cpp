
#include "reaction_kinematics.hpp"
#include "amplitude_sum.hpp"
#include "DsDpi.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void DsDpi_widths()
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
    DsDpi_swave swave(&kDsDpi, &Y, "S-wave");

    // D-wave recieves contributions from two diagrams: 
    // the triangle
    DsDpi_triangle triangle(&kDsDpi, &Y, "Triangle");

    // and the D1 tree transition
    DsDpi_tree tree(&kDsDpi, &Y, "Tree");

    // These can be summed up like this:
    amplitude_sum dwave(&kDsDpi, &Y, "D-wave");
    dwave.add({&tree, &triangle});

    // Both s and d-waves can then be added together 
    amplitude_sum full(&kDsDpi, &Y, "Full Sum");
    full.add({&dwave, &swave});

    // ---------------------------------------------------------------------------
    // Plotter object and options

    // Amplitudes to plot
    vector<amplitude*> amps;
    // amps.push_back(&full);
    amps.push_back(&swave);
    // amps.push_back(&dwave);
    // amps.push_back(&tree);
    // amps.push_back(&triangle);

    jpacGraph1D * plotter = new jpacGraph1D(0);

    int N = 100;
    bool PRINT = true;
    string xlabel, ylabel, filename;
    double xmin, xmax, ymin, ymax;
    ymin = 0.;

    // ---------------------------------------------------
    // AB channel 

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_b());
    xmax = (W - kDsDpi.mass_c());

    ymax = 7.;  

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}  [a.u.]";
    filename = "partial_width_ab.pdf";

    for (auto amp : amps)
    {
        auto dgammaab = [&](double w)
        {
            double sigma = w * w;
            return amp->dGamma_ab(s, sigma);
        };
        plotter->AddEntry(N, dgammaab, {xmin, xmax}, amp->get_id(), PRINT);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(0.7, 0.7);
    plotter->Plot(filename);

    // ---------------------------------------------------
    // BC channel 

    plotter->ClearData();

    xmin = (kDsDpi.mass_b() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_a());

    ymax = 7.;

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_b() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_b() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "partial_width_bc.pdf";

    for (auto amp : amps)
    {
        auto dgammabc = [&](double w)
        {
            double sigma = w * w;
            return amp->dGamma_bc(s, sigma);
        };
        plotter->AddEntry(N, dgammabc, {xmin, xmax}, amp->get_id(), PRINT);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(0.3, 0.7);
    plotter->Plot(filename);

    // ---------------------------------------------------
    // AC channel 

    plotter->ClearData();

    xmin = (kDsDpi.mass_a() + kDsDpi.mass_c());
    xmax = (W - kDsDpi.mass_b());

    ymax = 7.;

    xlabel   = "#sqrt{#sigma_{" + kDsDpi.particle_a() + kDsDpi.particle_c() + "}}   [GeV]";
    ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_c() + "}  [a.u.]";
    filename = "partial_width_ac.pdf";

    for (auto amp : amps)
    {
        auto dgammaac = [&](double w)
        {
            double sigma = w * w;
            return amp->dGamma_ac(s, sigma);
        };
        plotter->AddEntry(N, dgammaac, {xmin, xmax}, amp->get_id(), PRINT);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel/* , ymin, ymax */);
    plotter->SetLegend(0.3, 0.7);
    plotter->Plot(filename);

    // ---------------------------------------------------
    // Full width 

    plotter->ClearData();

    xmin = (4.1);
    xmax = (4.4);

    ymax = 1.;

    for (auto amp : amps)
    {
        auto gamma = [&](double w)
        {
            double s = w * w;
            return amp->Gamma(s);
        };
        plotter->AddEntry(N, gamma, {xmin, xmax}, amp->get_id(), PRINT);
    }

    xlabel   = "#sqrt{s}   [GeV]";
    ylabel   = "#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + kDsDpi.particle_c() + "}     [nb]";
    filename = "full_width.pdf";

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.3,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    delete plotter;
};  