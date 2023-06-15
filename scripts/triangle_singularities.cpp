#include "triangle.hpp"
#include "plotter.hpp"

void triangle_singularities()
{
    using namespace hadMolee;

    std::array<double,3> internal = {M_DSTAR, M_D1,     M_D};
    std::array<double,3> external = {M_Y4260, M_ZC3900, M_PION};

    triangle r_tri(triangle::kRelativistic, external, internal);
    r_tri.set_ieps(1E-3);

    triangle nr_tri(triangle::kNonrelativistic, external, internal);
    nr_tri.set_ieps(1E-3);

    std::vector<triangle> triangles = {r_tri, nr_tri};

    // ------------------------------------------------------
    
    // Plotter object
    plotter plotter;

    int N = 50;

    double xmin, xmax, ymin, ymax;
    string xlabel;

    std::string ylabel   = "#it{T}_{#it{D}*#it{D}_{1}#it{D}} #times 100";

    int i;
    auto re = [&](double E)
    {
        external[0] = E;
        triangles[i].set_external_masses(external);
        return real(triangles[i].eval() * 1.E2);
    };
    auto im = [&](double E)
    {
        external[0] = E;
        triangles[i].set_external_masses(external);
        return imag(triangles[i].eval() * 1.E2);
    };

    xmin = 4.2;
    xmax = 4.4;
    ymin = -1.;
    ymax = 5.;
    xlabel   = "#sqrt{#it{s}}  [GeV]";

    plot p1 = plotter.new_plot();
    p1.set_curve_points(N);
    i = 0;
    p1.add_curve({xmin, xmax}, re, "Real");
    i = 1;
    p1.add_dashed({xmin, xmax}, re);
    i = 0;
    p1.add_curve({xmin, xmax}, im, "Imaginary");
    i = 1;
    p1.add_dashed({xmin, xmax}, im);

    p1.add_header("#sqrt{#sigma_{D*D}} = 3.89 GeV");
    p1.set_ranges({xmin, xmax}, {ymin, ymax});
    p1.set_labels(xlabel, ylabel);
    p1.set_legend(0.7, 0.7);
    p1.save("tri_s.pdf");

    // Reset the masses
    external = {M_Y4260, M_ZC3900, M_PION};
    auto re2 = [&](double E)
    {
        external[1] = E;
        triangles[i].set_external_masses(external);
        return real(triangles[i].eval() * 1.E2);
    };
    auto im2 = [&](double E)
    {
        external[1] = E;
        triangles[i].set_external_masses(external);
        return imag(triangles[i].eval() * 1.E2);
    };

    xmin = 3.8;
    xmax = 4.1;
    ymin = 0.;
    ymax = 1.5;
    xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";

    plot p2 = plotter.new_plot();
    p2.set_curve_points(N);
    i = 0;
    p2.add_curve({xmin, xmax}, re2, "Real");
    i = 1;
    p2.add_dashed({xmin, xmax}, re2);
    i = 0;
    p2.add_curve({xmin, xmax}, im2, "Imaginary");
    i = 1;
    p2.add_dashed({xmin, xmax}, im2);
    
    p2.add_header("#sqrt{s} = 4.22 GeV");
    p2.set_ranges({xmin, xmax}, {ymin, ymax});
    p2.set_labels(xlabel, ylabel);
    p2.set_legend(0.7, 0.7);

    p2.save("tri_sig.pdf");
};  