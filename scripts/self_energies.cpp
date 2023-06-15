#include "molecule.hpp"
#include "Z_meson.hpp"
#include "Y_meson.hpp"
#include "breit_wigner.hpp"
#include "plotter.hpp"

void self_energies()
{   
    using namespace hadMolee; 
    using complex = std::complex<double>;

    // Dressed molecular propagator
    molecule Zc = make_molecule<DsD_molecule>();
    molecule Y  = make_molecule<D1D_molecule>();

    std::vector<molecule> mols = {Zc, Y};
    int i = 0;

    // ------------------------------------------------------

    double xmin, xmax, ymin, ymax;
    string xlabel, ylabel;

    // Plotter object
    plotter plotter;

    auto reSigma = [&] (double E)
    {
        return 100*real(mols[i]->self_energy(E*E));
    };
    auto imSigma = [&] (double E)
    {
        return 100*imag(mols[i]->self_energy(E*E));
    };

    auto G = [&] (double s)
    {
        std::array<double ,2> masses = mols[i]->constituent_masses();
        double m1 = masses[0];
        double m2 = masses[1];

        complex rho, xi;
        complex result;

        rho    = csqrt(Kallen(s, m1*m1, m2*m2)) / s;
        xi     = 1 - (m1+m2)*(m1+m2)/s;
        result = (rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1)) / PI;
        return result / (16*PI);
    };

    auto reG = [&] (double E)
    {
        return 100*real(G(E*E));
    };
    auto imG = [&] (double E)
    {
        return 100*imag(G(E*E));
    };

    xmin = 3.5;
    xmax = 4.1;
    ymin = -0.8;
    ymax = 2.3;

    xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";
    ylabel   = "#Sigma_{#it{D}*#it{D}} #times 100";

    i = 0;
    plot p1 = plotter.new_plot();
    p1.add_vertical(M_DSTAR + M_D);
    p1.add_curve({xmin, xmax}, reSigma, "Real");
    p1.add_dashed({xmin, xmax}, reG);
    p1.add_curve({xmin, xmax}, imSigma, "Imaginary");
    p1.add_dashed({xmin, xmax}, imG);

    p1.set_labels(xlabel, ylabel);
    p1.set_ranges({xmin, xmax}, {ymin, ymax});
    p1.set_legend( 0.25, 0.6);
    p1.save("z_sig.pdf");


    xmin = 4.0;
    xmax = 4.6;
    ymin = -0.8;
    ymax = 2.3;

    xlabel   = "#sqrt{#sigma_{#it{D}_{1}#it{D}}}  [GeV]";
    ylabel   = "#Sigma_{#it{D}_{1}*#it{D}} #times 100";


    i = 1;
    plot p2 = plotter.new_plot();
    p2.color_offset(2);
    p2.add_vertical(M_D1+ M_D);
    p2.add_curve({xmin, xmax}, reSigma, "Real");
    p2.add_dashed({xmin, xmax}, reG);
    p2.add_curve({xmin, xmax}, imSigma, "Imaginary");
    p2.add_dashed({xmin, xmax}, imG);

    p2.set_labels(xlabel, ylabel);
    p2.set_ranges({xmin, xmax}, {ymin, ymax});
    p2.set_legend( 0.25, 0.6);
    p2.save("y_sig.pdf");

    plotter.combine({2,1}, {p1,p2}, "self_energies.pdf");
};  