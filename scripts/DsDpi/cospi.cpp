#include "kinematics.hpp"
#include "plotter.hpp"
#include "DsDpi/D1.hpp"
#include "DsDpi/contact.hpp"
#include "Y(4260).hpp"

using namespace hadMolee;

void cospi()
{
    //--------------------------------------------------------------------------------------
    // Set up our kinematics

    kinematics kDsDpi      = make_kinematics({M_DSTAR, M_D, M_PION});

    //--------------------------------------------------------------------------------------
    // Set up our two linehsape models for the Y and psi(4160)

    lineshape  Y   = make_lineshape<Y_meson>();
    Y->set_parameters({4.2281, 42.6E-3, 3.115/sqrt(M_D1*M_D*4.23), 1./0.03946});

    lineshape  psi = make_lineshape<charmonium>();
    psi->set_parameters({4.191, 70.E-3, 1.});

    //--------------------------------------------------------------------------------------
    // Contact interactions

    amplitude  y_contact   = make_amplitude<DsDpi::contact>(kDsDpi, Y,   "Y contact");
    y_contact->set_parameters({21.733, -15.695});

    amplitude  psi_contact = make_amplitude<DsDpi::contact>(kDsDpi, psi, "#psi contact");
    psi_contact->set_parameters({0.627*4.66, -15.085});

    //--------------------------------------------------------------------------------------
    // Interactions involving intermediate D1's

    amplitude  tree        = make_amplitude<DsDpi::tree>       (kDsDpi, Y,   "D1 Tree");
    amplitude  triangle    = make_amplitude<DsDpi::triangle>   (kDsDpi, Y,   "D1 Triangle");

    //--------------------------------------------------------------------------------------
    // Add everything together
    amplitude full = y_contact + psi_contact + tree + triangle;

    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------

    auto plot_amp = [] (plot &p, amplitude amp, int option = hadMolee::kDefault, std::string id = "")
    {
        std::array<double,2> bounds = {-1., 1.};
        std::string label = (id == "") ? amp->get_id() : id;

        print("Plotting: ", label);
        divider(2);
        amp->set_option(option);
        double w = 4.23;
        auto xsection = [&] (double cos)
        {
            double x = amp->integrated_xsection(w*w, cos) * 1E3; // in pb
            print(cos, x);
            return x;
        };
        p.add_curve(bounds, xsection, label);
        line();
    };
    
    plotter plotter;
    plot sig = plotter.new_plot();
    sig.set_curve_points(50);
    sig.set_legend(0.15, 0.6);
    sig.set_labels("cos#theta_{#pi}", "#sigma [pb]");

    plot_amp(sig, full, 0, "Full Sum");

    // Save to file
    sig.save("cospi.pdf");
};