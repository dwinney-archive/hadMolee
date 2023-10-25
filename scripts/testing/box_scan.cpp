
#include "print.hpp"
#include "plotter.hpp"
#include "constants.hpp"
#include "box.hpp"
#include <memory>

void box_scan()
{
    using namespace hadMolee;
    using hadMolee::complex;

    // All on-shell meson masses are fixed and real for now
    double M_Y = 4.260;
    std::array<double, 4> internal = {M_D1, M_D, M_DSTAR, M_DSTAR}; 
    std::array<double, 4> external = {M_Y, M_PION, M_JPSI, M_PION};
    
    box TOPT(box::kTOPT);
    TOPT.set_external_masses(external);
    TOPT.set_internal_masses(internal);
    TOPT.set_max_calls(1E9);
    TOPT.add_width(0, W_D1);

    box LT(box::kLoopTools);
    LT.set_external_masses(external);
    LT.set_internal_masses(internal);
    LT.add_width(0, W_D1);

    // Fix pipi mass and scan the jpsi pi mass
    double s_pipi = pow(0.6, 2);
    std::array<double, 2> dalitz_bounds = {M_JPSI + M_PION + EPS,  M_Y - M_PION - EPS};
    auto reB = [&](box& implementation)
    {
        auto re = [&](double sqrts_psipi)
        {
            double s_psipi = sqrts_psipi * sqrts_psipi;
            double t = 2.*M_PION*M_PION + M_JPSI*M_JPSI + M_Y*M_Y - s_psipi - s_pipi;
            implementation.set_invariant_masses(s_psipi, t);
            return 100. * real(implementation.eval());
        };
        return re;
    };
    auto imB = [&](box& implementation)
    {
        return  [&](double sqrts_psipi)
        {
            double s_psipi = sqrts_psipi * sqrts_psipi;
            double t = 2.*M_PION*M_PION + M_JPSI*M_JPSI + M_Y*M_Y - s_psipi - s_pipi;
            implementation.set_invariant_masses(s_psipi, t);
            return 100. * imag(implementation.eval());
        };
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_curve_points(100);

    // p.set_ranges({3, 5}, {0., 0.6});
    p.set_legend(0.7, 0.7);
    p.add_header("#sqrt{#sigma_{#pi#pi}} = 0.6 GeV");
    p.set_labels("#sqrt{#sigma_{J/#psi#pi}}  [GeV]", "100 #times B_{D_{1}D*D*D}");
    p.add_vertical(M_JPSI + M_PION);
    p.add_vertical(M_Y - M_PION);

    p.add_curve( dalitz_bounds, reB(TOPT), "Real");
    p.add_dashed(dalitz_bounds, reB(LT));
    p.add_curve( dalitz_bounds, imB(TOPT), "Imaginary");
    p.add_dashed(dalitz_bounds, imB(LT));

    p.save("box.pdf");
};
