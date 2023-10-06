#include "molecule.hpp"
#include "lineshapes/Y(4260).hpp"
#include "breit_wigner.hpp"
#include "plotter.hpp"

void Y_propagator()
{   
    using namespace hadMolee; 
    using complex = std::complex<double>;

    // Dressed molecular propagator
    molecule Y = make_molecule<D1D_molecule>();

    // Nonrelativistic BW for comparison
    breit_wigner Y_BW(breit_wigner::kNonrelativistic, M_Y4260, W_Y4260);
    
    // ------------------------------------------------------

    double xmin = 4.0;
    double xmax = 4.4;

    double ymin = -20.;
    double ymax =  40.;

    string xlabel   = "#it{E} [GeV]";
    string ylabel   = "#it{i}#it{G}_{#it{Y}}";
    string filename = "Z_prop.pdf";

    // Plotter object
    plotter plotter;

    auto reiG = [&] (double E)
    {
        return real(I*Y->propagator(E*E));
    };
    auto imiG = [&] (double E)
    {
        return imag(I*Y->propagator(E*E));
    };

    auto reiBW = [&] (double E)
    {
        return real(I*Y_BW.eval(E));
    };
    auto imiBW = [&] (double E)
    {
        return imag(I*Y_BW.eval(E));
    };

    plot p = plotter.new_plot();
    p.set_curve_points(500);
    p.add_curve( {xmin, xmax}, reiG, "Real");
    p.add_dashed({xmin, xmax}, reiBW);

    p.add_curve( {xmin, xmax}, imiG, "Imaginary");
    p.add_dashed({xmin, xmax}, imiBW);
    p.add_vertical(M_D1 + M_D);
    
    p.set_ranges({xmin, xmax}, {ymin, ymax});
    p.set_labels(xlabel, ylabel);
    p.set_legend(0.4, 0.75);
    p.save("Z_prop.pdf");
};  