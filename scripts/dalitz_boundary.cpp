
#include "reaction_kinematics.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void dalitz_boundary()
{   
    reaction_kinematics kDsDpi({M_DSTAR, M_D, M_PION});
    kDsDpi.set_particle_labels("D*", "D", "#pi");

    double W = 4.260;
    double s = W*W;

    int N = 200;
    bool PRINT = true;

    double xmin = (kDsDpi.mass_a() + kDsDpi.mass_b());
    double xmax = (W - kDsDpi.mass_c());

    double ymin = 2.;
    double ymax = 2.3;

    string label_ab = "#sqrt{#sigma_{"+ kDsDpi.particle_a() + kDsDpi.particle_b() + "}}   [GeV]";
    string label_bc = "#sqrt{#sigma_{"+ kDsDpi.particle_b() + kDsDpi.particle_c() + "}}   [GeV]";
    string label_ac = "#sqrt{#sigma_{"+ kDsDpi.particle_a() + kDsDpi.particle_c() + "}}   [GeV]";

    string filename = "dalitz_test.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto sbc_min = [&](double wab)
    {
        double sab = wab * wab;
        return sqrt(kDsDpi.sbc_from_sab(s, sab, -1.));
    };
    auto sbc_max = [&](double wab)
    {
        double sab = wab * wab;
        return sqrt(kDsDpi.sbc_from_sab(s, sab, +1));
    };

    plotter->AddEntry(N, sbc_max, {xmin, xmax}, "", PRINT);
    plotter->AddEntry(N, sbc_min, {xmin, xmax}, "", PRINT);
    plotter->SetXaxis(label_ab, xmin, xmax);
    plotter->SetYaxis(label_bc, ymin, ymax);
    plotter->SetLegend(false);
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;

};  