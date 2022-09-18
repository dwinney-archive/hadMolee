
#include "reaction_kinematics.hpp"
#include "amplitude.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void width_phase_space()
{   
    
    double W = 4.260;
    double s = W*W;

    reaction_kinematics kDsDpi({M_DSTAR, M_D, M_PION});
    kDsDpi.set_particle_labels({"#it{D}*", "#it{D}", "#pi"});

    phase_space DsDpi(&kDsDpi, "phase space");
    DsDpi.normalize(1., s);

    int N = 100;
    bool PRINT = true;

    double xmin = (kDsDpi.mass_a() + kDsDpi.mass_b());
    double xmax = (W - kDsDpi.mass_c());

    double ymin = 0.;
    double ymax = 1.;

    string label_ab = "#sqrt{#sigma_{"+ kDsDpi.particle_a() + kDsDpi.particle_b() + "}}   [GeV]";
    string ylabel   = "#it{d}#Gamma_{" + kDsDpi.particle_a() + kDsDpi.particle_b() + "}";
    string filename = "widthab.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto dgamma = [&](double wab)
    {
        double sab = wab * wab;
        return DsDpi.dGamma_ab(s, sab);
    };

    plotter->AddEntry(N, dgamma, {xmin, xmax}, "", PRINT);
    plotter->SetXaxis(label_ab, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(false);
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;

};  