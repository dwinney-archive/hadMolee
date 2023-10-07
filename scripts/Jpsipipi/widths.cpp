#include "kinematics.hpp"
#include "lineshapes/Y(4260).hpp"
#include "Jpsipipi/boxes.hpp"
#include "plotter.hpp"

using namespace hadMolee;

void widths()
{
    // Masses
    double ma = M_JPSI;
    double mb = M_PION;
    double mc = M_PION;

    kinematics kJpsipipi  = make_kinematics({ma, mb, mc});
    lineshape  Y          = make_lineshape<D1D_molecule>();

    amplitude  box_I      = make_amplitude<Jpsipipi::Box_I>(kJpsipipi, Y, "Box I");
    // amplitude box_I = make_amplitude<phase_space>(kJpsipipi, Y, "");

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
        divider(2);
        int i = 0;
        auto xsection = [&] (double E)
        {
            double x = amp->differential_xsection(abc, W*W, E*E) * 1E3; // in pb
            print(i, E, x);
            i++;
            return x;
        };
        p.add_curve(bounds, xsection, amp->get_id());
        line();
    };

    plotter plotter;

    plot pab = plotter.new_plot();
    pab.set_curve_points(50);
    pab.set_labels("#it{E}_{#it{J}/#psi#pi}   [GeV]", "d#sigma/dE^{2} [nb / GeV^{-2}]");
    plot_amp(pab, ab, box_I);
    pab.save("widths.pdf");

    // plot pbc = plotter.new_plot();
    // pbc.set_curve_points(10);
    // pbc.set_labels("#it{E}_{#pi#pi}   [GeV]", "d#sigma/dE^{2}  [nb / GeV^{-2}]");
    // plot_amp(pbc, bc, box_I);

    // plot pac = plotter.new_plot();
    // pac.set_curve_points(10);
    // pac.set_labels("#it{E}_{#it{J}/#psi#pi}  [GeV]", "d#sigma/dE^{2}  [nb / GeV^{-2}]");
    // plot_amp(pac, ac, box_I);

    // plotter.combine({1,3}, {pab, pac, pbc}, "widths.pdf");
};