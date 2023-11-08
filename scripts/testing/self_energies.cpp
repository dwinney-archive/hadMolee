#include "molecule.hpp"
#include "Y(4260).hpp"
#include "breit_wigner.hpp"
#include "plotter.hpp"

void self_energies()
{   
    using namespace hadMolee; 
    using complex = std::complex<double>;

    molecule Y  = make_molecule(M_D1, M_D);
    Y->set_parameters({4.2281, 42.6E-3, 0.1145/sqrt(4.2281*M_D1*M_DSTAR)});

    // ------------------------------------------------------

    double xmin, xmax, ymin, ymax;
    string xlabel, ylabel;

    // Plotter object
    plotter plotter;

    auto reSigma = [&] (double E)
    {
        return 100*real(Y->self_energy(E*E));
    };
    auto imSigma = [&] (double E)
    {
        return 100*imag(Y->self_energy(E*E));
    };

    auto G = [&] (double s)
    {
        std::array<double ,2> masses = Y->constituent_masses();
        double m1 = masses[0];
        double m2 = masses[1];

        if (sqrt(s) < m1 + m2) return 0.;
        return - sqrt(Kallen(s, m1*m1, m2*m2)) / s / (16*PI);
    };

    auto reG = [&] (double E)
    {
        return 100*real(G(E*E));
    };
    auto imG = [&] (double E)
    {
        return 100*imag(G(E*E));
    };

    xmin = 4.0;
    xmax = 4.6;
    ymin = -1.8;
    ymax = 2.3;

    xlabel   = "#sqrt{#sigma_{#it{D}_{1}#it{D}}}  [GeV]";
    ylabel   = "#Sigma_{#it{D}_{1}*#it{D}} #times 100";

    plot p2 = plotter.new_plot();
    p2.set_labels(xlabel, ylabel);
    p2.set_ranges({xmin, xmax}, {ymin, ymax});
    p2.set_legend( 0.25, 0.6);

    p2.add_vertical(M_D1+ M_D);
    p2.add_curve({xmin, xmax}, reSigma, "Real");
    p2.add_curve({xmin, xmax}, imSigma, "Imaginary");
    p2.color_offset(2);
    p2.add_dashed({xmin, xmax}, reG);

    p2.save("y_sigma.pdf");
};  