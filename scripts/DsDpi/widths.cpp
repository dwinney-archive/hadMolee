#include "kinematics.hpp"
#include "plotter.hpp"
#include "DsDpi.hpp"
#include "lineshapes/Y(4260).hpp"

using namespace hadMolee;

void widths()
{
    // Masses
    double ma = M_DSTAR;
    double mb = M_D;
    double mc = M_PION;

    kinematics kDsDpi    = make_kinematics({ma, mb, mc});
    lineshape  Y         = make_lineshape<D1D_molecule>();

    amplitude  swave     = make_amplitude<DsDpi::swave>(    kDsDpi, Y, "S-wave");
    amplitude  tree      = make_amplitude<DsDpi::tree>(     kDsDpi, Y, "Tree");
    amplitude  triangle  = make_amplitude<DsDpi::triangle>( kDsDpi, Y, "Triangle");

    amplitude  dwave     = tree + triangle;
    dwave->set_id("D-wave");

    amplitude  sum     = swave + dwave;
    sum->set_id("Sum");

    double W = 4.23;

    auto plot_amp = [W, ma, mb, mc] (plot &p, subchannel abc, amplitude amp)
    {
        std::array<double,2> bounds;
        switch (abc)
        {
            case ab: bounds = {ma + mb, W - mc}; break;
            case bc: bounds = {mc + mb, W - ma}; break;
            case ac: bounds = {ma + mc, W - mb}; break;
        }

        print("Plotting: ", amp->get_id());
        auto xsection = [&] (double E)
        {
            double x = amp->differential_xsection(abc, W*W, E*E) * 1E3; // in pb
            print(E, x);
            return x;
        };
        p.add_curve(bounds, xsection, amp->get_id());
        line();
    };

    plotter plotter;
    plot pab = plotter.new_plot();
    pab.set_curve_points(100);
    pab.set_labels("#it{E}_{#it{D}#it{D}*}   [GeV]", "d#sigma/dE^{2} [nb / GeV^{-2}]");

    plot_amp(pab, ab, swave);
    plot_amp(pab, ab, dwave);
    plot_amp(pab, ab, sum);

    plot pbc = plotter.new_plot();
    pbc.set_curve_points(100);
    pbc.set_labels("#it{E}_{#pi#it{D}}   [GeV]", "d#sigma/dE^{2}  [nb / GeV^{-2}]");

    plot_amp(pbc, bc, swave);
    plot_amp(pbc, bc, dwave);
    plot_amp(pbc, bc, sum);

    plot pac = plotter.new_plot();
    pac.set_curve_points(100);
    pac.set_labels("#it{E}_{#pi#it{D}*}   [GeV]", "d#sigma/dE^{2}  [nb / GeV^{-2}]");

    plot_amp(pac, ac, swave);
    plot_amp(pac, ac, dwave);
    plot_amp(pac, ac, sum);

    plotter.combine({1,3}, {pab, pac, pbc}, "widths.pdf");
};