#include "molecule.hpp"
#include "lineshapes/Z(3900).hpp"
#include "breit_wigner.hpp"
#include "plotter.hpp"

void Z_propagator()
{   
    using namespace hadMolee; 
    using complex = std::complex<double>;

    // Dressed molecular propagator
    molecule Zc = make_molecule<DsD_molecule>();

    // Nonrelativistic BW for comparison
    breit_wigner Zc_BW(breit_wigner::kNonrelativistic, M_ZC3900, W_ZC3900);
    
    // ------------------------------------------------------

    double xmin = 3.8;
    double xmax = 4.1;

    double ymin = -20.;
    double ymax =  40.;

    string xlabel   = "#it{E}  [GeV]";
    string ylabel   = "#it{i}#it{G}_{#it{Z}}";
    string filename = "Z_prop.pdf";

    // Plotter object
    plotter plotter;

    auto reiG = [&] (double E)
    {
        return real(I*Zc->propagator(E*E));
    };
    auto imiG = [&] (double E)
    {
        return imag(I*Zc->propagator(E*E));
    };

    auto reiBW = [&] (double E)
    {
        return real(I*Zc_BW.eval(E));
    };
    auto imiBW = [&] (double E)
    {
        return imag(I*Zc_BW.eval(E));
    };

    plot p = plotter.new_plot();
    p.color_offset(2);
    p.set_curve_points(500);
    p.add_curve( {xmin, xmax}, reiG, "Real");
    p.add_dashed({xmin, xmax}, reiBW);

    p.add_curve( {xmin, xmax}, imiG, "Imaginary");
    p.add_dashed({xmin, xmax}, imiBW);
    p.add_vertical(M_DSTAR + M_D);
    
    p.set_ranges({xmin, xmax}, {ymin, ymax});
    p.set_labels(xlabel, ylabel);
    p.set_legend(0.7, 0.75);
    p.save("Z_prop.pdf");
};  