
#include "print.hpp"
#include "plotter.hpp"
#include "constants.hpp"
#include "box.hpp"

void box_scan()
{
    using namespace hadMolee;
    using hadMolee::complex;

    // All on-shell meson masses are fixed and real for now
    double M_Y = 4.260;
    std::array<double, 4> internal = {M_D1, M_D, M_DSTAR, M_DSTAR}; 
    std::array<double, 4> external = {M_Y, M_PION, M_JPSI, M_PION};
    
    // Use the brute-force box first
    box Box(box::kLoopTools);
    Box.set_external_masses(external);
    Box.set_internal_masses(internal);

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_curve_points(200);

    // Fix pipi mass and scan the jpsi pi mass
    double s_pipi = pow(0.6, 2);

    auto B = [&](double sqrts_psipi)
    {
        double s_psipi = sqrts_psipi * sqrts_psipi;
        double t = 2.*M_PION*M_PION + M_JPSI*M_JPSI + M_Y*M_Y - s_psipi - s_pipi;
        Box.set_invariant_masses(s_psipi, t);
        return 100. * Box.eval();
    };
    
    auto reB = [&](double w){ return std::real(B(w)); };
    auto imB = [&](double w){ return std::imag(B(w)); };
    
    print("Thresholds");
    print("D1 D* = ", M_D1 + M_DSTAR);
    print("D* D  = ", M_DSTAR + M_D);
    print("D* D* = ", 2*M_DSTAR);
    print("D1 D  = ", M_D1 + M_D);
    print("Jpsi Pi pi =", M_JPSI + 2*M_PION);

    p.set_ranges({3, 5}, {0., 0.6});
    p.set_legend(0.7, 0.7);
    p.add_header("#sqrt{#sigma_{#pi#pi}} = 0.6 GeV");
    p.set_labels("#sqrt{#sigma_{J/#psi#pi}}  [GeV]", "100 #times B_{D_{1}D*D*D}");
    p.add_vertical(M_JPSI + M_PION);
    p.add_vertical(M_Y - M_PION);
    p.add_curve({3., 5.}, reB, "Real");
    Box.add_width(0, W_D1);
    p.add_dashed({3., 5.}, reB);
    Box.add_width(0, 0.);
    p.add_curve({3,  5.}, imB, "Imag");
    Box.add_width(0, W_D1);
    p.add_dashed({3., 5.}, imB);

    p.save("box.pdf");
};
