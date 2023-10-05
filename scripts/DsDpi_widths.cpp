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
        return to_plot->differential_xsection(ab, s, Eab*Eab);
    };
    auto dGamma_ac = [&] (double Eac)
    {
        return to_plot->differential_xsection(ac, s, Eac*Eac);
    };
    auto dGamma_bc = [&] (double Ebc)
    {
        return to_plot->differential_xsection(bc, s, Ebc*Ebc);
    };

    plotter plotter;
    plot pab = plotter.new_plot();
    pab.set_curve_points(100);
    pab.set_labels("#it{E}_{#it{D}#it{D}*}   [GeV]", "d#sig [nb / GeV^{-2}]");

    to_plot = swave;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());
    to_plot = dwave;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());
    to_plot = sum;
    pab.add_curve({ma + mb, W - mc}, dGamma_ab, to_plot->get_id());

    plot pbc = plotter.new_plot();
    pbc.set_curve_points(100);
    pbc.set_labels("#it{E}_{#pi#it{D}}   [GeV]", "d#sigma [nb / GeV^{-2}]");
    to_plot = swave;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());
    to_plot = dwave;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());
    to_plot = sum;
    pbc.add_curve({mc + mb, W - ma}, dGamma_bc, to_plot->get_id());

    plot pac = plotter.new_plot();
    pac.set_curve_points(100);
    pac.set_labels("#it{E}_{#pi#it{D}*}   [GeV]", "d#sigma [nb / GeV^{-2}]");
    to_plot = swave;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());
    to_plot = dwave;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());
    to_plot = sum;
    pac.add_curve({ma + mc, W - mb}, dGamma_ac, to_plot->get_id());


    plotter.combine({3,1}, {pab, pac, pbc}, "Gammas.pdf");
};