#include "molecule.hpp"
#include "Y(4260).hpp"
#include "breit_wigner.hpp"
#include "plotter.hpp"

void Y_propagator()
{   
    using namespace hadMolee; 
    using complex = std::complex<double>;

    // Dressed molecular propagator
    lineshape Y = make_lineshape<Y_meson>();
    Y->set_parameters({4.2281, 42.6E-3, 0.1145/sqrt(4.2281*M_D1*M_DSTAR), 0.672});

    // Nonrelativistic BW for comparison
    breit_wigner Y_BW(breit_wigner::kNonrelativistic, M_Y4260, W_Y4260);
    
    // ------------------------------------------------------

    double xmin = 4.0;
    double xmax = 4.4;

    double ymin = -10.;
    double ymax =  4.;

    string xlabel   = "#it{E} [GeV]";
    string ylabel   = "#it{G}_{#it{Y}}";
    string filename = "Z_prop.pdf";

    // Plotter object
    plotter plotter;

    auto reiG = [&] (double E)
    {
        print(E, real(Y->propagator(E*E)));
        return real(Y->propagator(E*E));
    };
    auto imiG = [&] (double E)
    {
        print(E, imag(Y->propagator(E*E)));
        return imag(Y->propagator(E*E));
    };

    plot p = plotter.new_plot();
    p.set_curve_points(100);
    p.add_curve( {xmin, xmax}, reiG, "Real");
    
    line();
    line();
    p.add_curve( {xmin, xmax}, imiG, "Imaginary");
    p.add_vertical(M_D1 + M_D);
    
    p.set_ranges({xmin, xmax}, {ymin, ymax});
    p.set_labels(xlabel, ylabel);
    p.set_legend(0.4, 0.75);
    p.save("Z_prop.pdf");
};  