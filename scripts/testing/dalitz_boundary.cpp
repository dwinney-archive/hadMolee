
#include "kinematics.hpp"
#include "plotter.hpp"

void dalitz_boundary()
{   
    using namespace hadMolee;

    kinematics kDsDpi = make_kinematics({M_DSTAR, M_D, M_PION});

    double W = 4.260;
    double s = W*W;

    double xmin = (kDsDpi->mass_a() + kDsDpi->mass_b());
    double xmax = (W - kDsDpi->mass_c());

    double ymin = 2.;
    double ymax = 2.3;

    string label_ab = "#sqrt{#sigma_{"+ kDsDpi->particle_a() + kDsDpi->particle_b() + "}}   [GeV]";
    string label_bc = "#sqrt{#sigma_{"+ kDsDpi->particle_b() + kDsDpi->particle_c() + "}}   [GeV]";
    string label_ac = "#sqrt{#sigma_{"+ kDsDpi->particle_a() + kDsDpi->particle_c() + "}}   [GeV]";

    string filename = "dalitz_test.pdf";

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto sbc_min = [&](double wab)
    {
        double sab = wab * wab;
        return sqrt(kDsDpi->sbc_from_sab(s, sab, -1));
    };
    auto sbc_max = [&](double wab)
    {
        double sab = wab * wab;
        return sqrt(kDsDpi->sbc_from_sab(s, sab, +1));
    };

    plotter plotter; 
    
    plot p = plotter.new_plot();
    p.set_curve_points(200);
    p.set_legend(false);
    p.set_ranges({xmin, xmax}, {ymin, ymax});
    p.set_labels(label_ab, label_bc);

    p.add_curve({xmin, xmax}, sbc_max);
    p.add_curve({xmin, xmax}, sbc_min);

    p.save(filename);
};  