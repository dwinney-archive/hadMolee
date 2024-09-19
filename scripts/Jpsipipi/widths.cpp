#include "kinematics.hpp"
#include "plotter.hpp"
#include "Y(4260).hpp"
#include "Jpsipipi/boxes.hpp"
#include "Jpsipipi/triangles.hpp"
// #include "Jpsipipi/counter_terms.hpp"


void widths()
{
    using namespace hadMolee;
    
    // Masses
    double ma = M_JPSI;
    double mb = M_PION;
    double mc = M_PION;

    kinematics kJpsipipi  = make_kinematics({ma, mb, mc});
    lineshape  Y          = make_lineshape<Y_meson>();
    Y->set_parameters({4.2281, 42.6E-3, 3.115/sqrt(M_D1*M_D0*4.2281), 1./0.03946});

    amplitude  box_I      = make_amplitude<Jpsipipi::box_I>(  kJpsipipi, Y, "Box I");
    amplitude  box_II     = make_amplitude<Jpsipipi::box_II>( kJpsipipi, Y, "Box II");
    amplitude  box_III    = make_amplitude<Jpsipipi::box_III>(kJpsipipi, Y, "Box III");
    // amplitude  box_sum = box_I + box_II + box_III;

    amplitude tri_I       = make_amplitude<Jpsipipi::triangle_I>(  kJpsipipi, Y, "Triangle I");
    amplitude tri_II      = make_amplitude<Jpsipipi::triangle_II>( kJpsipipi, Y, "Triangle II");
    amplitude tri_III     = make_amplitude<Jpsipipi::triangle_III>(kJpsipipi, Y, "Triangle III");
    amplitude tri_sum     = tri_I + tri_II + tri_III;

    // amplitude tri_CT     = make_amplitude<Jpsipipi::triangle_CT>(kJpsipipi, Y, "Triangle CT");

    double W = 4.23;
    auto plot_amp = [&] (plot &p, subchannel abc, amplitude amp)
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
        auto xsection = [&] (double sig)
        {
            double x = 2*W*amp->differential_xsection(abc, W*W, sig*sig) * 1E3; // in pb
            print(sig, x);
            return x;
        };
        p.add_curve(bounds, xsection);
        line();
    };

    plotter plotter;

    double E   = 4.23;
    double Eab = 3.4;
    double Ebc = 0.8;
    double cos = 1.0;


    // print(tri_III->differential_xsection(E*E, Eab*Eab, Ebc*Ebc, cos));
    // print(box_I->differential_xsection(E*E, Eab*Eab, Ebc*Ebc, cos));

    // Jpsi pi system
    plot pab = plotter.new_plot();
    pab.set_curve_points(50);
    pab.set_labels("#it{E}_{#it{J}/#psi#pi}   [GeV]", """d#sigma/dE [pb / GeV^{-1}]");
    // plot_amp(pab, ab, tri_II);
    // plot_amp(pab, ab, tri_I);
    // plot_amp(pab, ab, tri_III);
    // plot_amp(pab, ab, tri_I + tri_II + tri_III);
    // plot_amp(pab, ab, box_II);

    // box_I->set_option(1);
    // plot_amp(pab, ab, box_I);
    // box_I->set_option(2);
    // plot_amp(pab, ab, box_I);
    box_I->set_option(0);
    plot_amp(pab, ab, box_I);

    // plot_amp(pab, ab, box_III);
    pab.save("widths.pdf");
};