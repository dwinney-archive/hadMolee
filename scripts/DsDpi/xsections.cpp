#include "kinematics.hpp"
#include "plotter.hpp"
#include "DsDpi.hpp"
#include "lineshapes/Y(4260).hpp"
#include "lineshapes/charmonium.hpp"

using namespace hadMolee;

void xsections()
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

    lineshape psi4160    = make_lineshape<charmonium>({4.191, 70.E-3, 1.});

    amplitude  dwave     = tree + triangle;
    dwave->set_id("D-wave");

    amplitude  sum     = swave + dwave;
    sum->set_id("Sum");

    amplitude to_plot;

    // Lambdas for each of the partial widths
    auto xsection = [&] (double w)
    {
        double x = to_plot->integrated_xsection(w*w) * 1E3;
        print(w, x);
        return x;
    };
    std::array<double,2> bounds = {4.05, 4.40};

    plotter plotter;
    plot sig = plotter.new_plot();
    sig.set_curve_points(15);
    sig.set_labels("#sqrt{#it{s}}  [GeV]", "#sigma [pb]");

    to_plot = swave;
    sig.add_curve(bounds, xsection, to_plot->get_id());
    to_plot = dwave;
    sig.add_curve(bounds, xsection, to_plot->get_id());
    to_plot = sum;
    sig.add_curve(bounds, xsection, to_plot->get_id());

    sig.save("sig.pdf");
};