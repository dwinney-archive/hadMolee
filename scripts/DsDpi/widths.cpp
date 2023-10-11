#include "kinematics.hpp"
#include "lineshape.hpp"
#include "plotter.hpp"
#include "DsDpi/D1.hpp"
#include "DsDpi/contact.hpp"
#include "Y(4260).hpp"

using namespace hadMolee;

void widths()
{
    //--------------------------------------------------------------------------------------
    // Set up our kinematics

    kinematics kDsDpi      = make_kinematics({M_DSTAR, M_D, M_PION});

    //--------------------------------------------------------------------------------------
    // Set up our two linehsape models for the Y and psi(4160)

    lineshape  Y   = make_lineshape<Y_meson>();
    Y->set_parameters({4.2281, 42.6E-3, 0.1145/sqrt(4.2281*M_D1*M_DSTAR), 0.672});

    lineshape  psi = make_lineshape<charmonium>();
    psi->set_parameters({4.191, 70.E-3, 6.9E-6});

    //--------------------------------------------------------------------------------------
    // Contact interactions

    amplitude  y_contact   = make_amplitude<DsDpi::contact>(kDsDpi, Y,   "Y contact");
    y_contact->set_parameters({1.110, -15.5});

    amplitude  psi_contact = make_amplitude<DsDpi::contact>(kDsDpi, psi, "#psi contact");
    psi_contact->set_parameters({0.7/psi->photon_coupling(), -14.96});

    //--------------------------------------------------------------------------------------
    // Interactions involving intermediate D1's

    amplitude  tree        = make_amplitude<DsDpi::tree>       (kDsDpi, Y,   "D1 Tree");
    amplitude  triangle    = make_amplitude<DsDpi::triangle>   (kDsDpi, Y,   "D1 Triangle");

    //--------------------------------------------------------------------------------------
    // Plot results

    // Fixed e+ e- energy
    double W = 4.23;

    auto plot_amp = [W, kDsDpi] (plot &p, subchannel abc, amplitude amp, int option = hadMolee::kDefault, std::string id = "")
    {
        std::array<double,2> bounds;
        switch (abc)
        {
            case ab: bounds = {kDsDpi->mass_a() + kDsDpi->mass_b(), W - kDsDpi->mass_c()}; break;
            case bc: bounds = {kDsDpi->mass_b() + kDsDpi->mass_c(), W - kDsDpi->mass_a()}; break;
            case ac: bounds = {kDsDpi->mass_a() + kDsDpi->mass_c(), W - kDsDpi->mass_b()}; break;
        }
        std::string label = (id == "") ? amp->get_id() : id;
        amp->set_option(option);

        print("Plotting: ", label);
        divider(2);
        auto xsection = [&] (double E)
        {
            double x = amp->differential_xsection(abc, W*W, E*E); // in pb
            print(E, x);
            return x;
        };
        p.add_curve(bounds, xsection, label);
        line();
    };

    plotter plotter;
    plot pab = plotter.new_plot();
    pab.set_curve_points(100);
    pab.set_legend(0.8, 0.7);
    pab.set_labels("#it{E}_{#it{D}#it{D}*}   [GeV]", "d#sigma/dE^{2} [nb / GeV^{-2}]");
    // plot_amp(pab, ab, y_contact);
    plot_amp(pab, ab, psi_contact);
    plot_amp(pab, ab, tree,     DsDpi::D1::kSwaveOnly, "S-wave Tree");
    plot_amp(pab, ab, tree,     DsDpi::D1::kDwaveOnly, "D-wave Tree");
    plot_amp(pab, ab, triangle, DsDpi::D1::kSwaveOnly, "S-wave Triangle");
    plot_amp(pab, ab, triangle, DsDpi::D1::kDwaveOnly, "D-Wave Triangle");

    plot pbc = plotter.new_plot();
    pbc.set_curve_points(100);
    pbc.set_legend(0.8, 0.7);
    pbc.set_labels("#it{E}_{#it{D}#pi}   [GeV]", "d#sigma/dE^{2} [nb / GeV^{-2}]");
    // plot_amp(pbc, bc, y_contact);
    plot_amp(pbc, bc, psi_contact);
    plot_amp(pbc, bc, tree,     DsDpi::D1::kSwaveOnly, "S-wave Tree");
    plot_amp(pbc, bc, tree,     DsDpi::D1::kDwaveOnly, "D-wave Tree");
    plot_amp(pbc, bc, triangle, DsDpi::D1::kSwaveOnly, "S-wave Triangle");
    plot_amp(pbc, bc, triangle, DsDpi::D1::kDwaveOnly, "D-Wave Triangle");

    plot pac = plotter.new_plot();
    pac.set_curve_points(100);
    pac.set_legend(0.8, 0.7);
    pac.set_labels("#it{E}_{#it{D}*#pi}   [GeV]", "d#sigma/dE^{2} [nb / GeV^{-2}]");
    // plot_amp(pac, ac, y_contact);
    plot_amp(pac, ac, psi_contact);
    plot_amp(pac, ac, tree,     DsDpi::D1::kSwaveOnly, "S-wave Tree");
    plot_amp(pac, ac, tree,     DsDpi::D1::kDwaveOnly, "D-wave Tree");
    plot_amp(pac, ac, triangle, DsDpi::D1::kSwaveOnly, "S-wave Triangle");
    plot_amp(pac, ac, triangle, DsDpi::D1::kDwaveOnly, "D-Wave Triangle");

    // Save to file
    plotter.combine({1,3}, {pab, pbc, pac}, "widths.pdf");
};