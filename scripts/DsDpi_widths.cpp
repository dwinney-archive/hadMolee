#include "kinematics.hpp"
#include "DsDpi.hpp"
#include "Y_meson.hpp"
#include "plotter.hpp"

using namespace hadMolee;

void DsDpi_widths()
{
    // Masses
    double ma = M_DSTAR;
    double mb = M_D;
    double mc = M_PION;

    kinematics kDsDpi    = make_kinematics({ma, mb, mc});
    lineshape  Y         = make_lineshape<D1D_molecule>();

    amplitude  swave     = make_amplitude<DsDpi_swave>(    kDsDpi, Y, "S-wave");
    amplitude  tree      = make_amplitude<DsDpi_tree>(     kDsDpi, Y, "Tree");
    amplitude  triangle  = make_amplitude<DsDpi_triangle>( kDsDpi, Y, "Triangle");

    amplitude  dwave     = tree + triangle;
    dwave->set_id("D-wave");

    amplitude  sum     = swave + dwave;
    sum->set_id("Sum");

    amplitude to_plot;

    double W = 4.23;
    double s = W*W;

    // Lambdas for each of the partial widths
    auto dGamma_ab = [&] (double Eab)
    {
        return to_plot->dGamma(ab, s, Eab*Eab) * 1E3;
    };
    auto dGamma_ac = [&] (double Eac)
    {
        return to_plot->dGamma(ac, s, Eac*Eac) * 1E3;
    };
    auto dGamma_bc = [&] (double Ebc)
    {
        return to_plot->dGamma(bc, s, Ebc*Ebc) * 1E3;
    };

    plotter plotter;
    plot pab = plotter.new_plot();
    pab.set_curve_points(100);
    pab.set_labels("#it{E}_{#it{D}^{*}#it{D}}   [GeV]", "#it{d}#Gamma #times 10^{3}   [MeV]");

    to_plot = dwave;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());
    to_plot = triangle;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());
    to_plot = tree;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());

    plot pbc = plotter.new_plot();
    pbc.set_curve_points(100);
    pbc.set_labels("#it{E}_{#it{D}#pi}   [GeV]", "#it{d}#Gamma #times 10^{3}   [MeV]");
    to_plot = dwave;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());
    to_plot = triangle;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());
    to_plot = tree;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());

    plot pac = plotter.new_plot();
    pac.set_curve_points(100);
    pac.set_labels("#it{E}_{#it{D}^{*}#pi}   [GeV]", "#it{d}#Gamma #times 10^{3}   [MeV]");
    to_plot = dwave;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());
    to_plot = triangle;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());
    to_plot = tree;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());


    plotter.combine({3,1}, {pab, pac, pbc}, "Gammas.pdf");
};
