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

    lineshape  psi4160   = make_lineshape<charmonium>({4.191, 70.E-3, 6.9E-6}, "#psi");
    amplitude  psi       = make_amplitude<DsDpi::psi_contact>(kDsDpi, psi4160, "#psi(4160)");

    amplitude  dwave     = tree + triangle;
    dwave->set_id("D-wave");

    amplitude  sum       = swave + dwave;
    sum->set_id("Sum");

    auto plot_amp = [] (plot &p, amplitude amp)
    {
        std::array<double,2> bounds = {4.05, 4.40};

        print("Plotting: ", amp->get_id());
        auto xsection = [&] (double w)
        {
            double x = amp->integrated_xsection(w*w) * 1E3; // in pb
            print(w, x);
            return x;
        };
        p.add_curve(bounds, xsection, amp->get_id());
        line();
    };


    plotter plotter;
    plot sig = plotter.new_plot();
    sig.set_curve_points(50);
    sig.set_labels("#sqrt{#it{s}}  [GeV]", "#sigma [pb]");

    // plot_amp(sig, swave);
    plot_amp(sig, psi);
    // plot_amp(sig, dwave);
    // plot_amp(sig, sum);

    sig.save("xsections.pdf");
};