#include "triangle.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void triangle_singularities()
{   
    array<double,3> internal = {M_DSTAR, M_D1,     M_D};
    array<double,3> external = {M_Y4260, M_ZC3900, M_PION};

    relativistic_triangle r_triangle;
    r_triangle.set_internal_masses(internal);
    r_triangle.set_external_masses(external);
    r_triangle.set_ieps(1.E-5);

    nonrelativistic_triangle nr_triangle;
    nr_triangle.set_internal_masses(internal);
    nr_triangle.set_external_masses(external);
    nr_triangle.set_ieps(1.E-5);

    vector<triangle*> triangles = {&r_triangle, &nr_triangle};
    
    // ------------------------------------------------------
    int N = 100;
    bool PRINT = true;

    double xmin, xmax, ymin, ymax;
    string xlabel, filename;

    string ylabel   = "#it{T}_{#it{D}*#it{D}_{1}#it{D}} #times 100";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D(0);

    int i;
    auto re = [&](double E)
    {
        external[0] = E;
        triangles[i]->set_external_masses(external);
        return real(triangles[i]->eval() * 1.E2);
    };
    auto im = [&](double E)
    {
        external[1] = E;
        triangles[i]->set_external_masses(external);
        return imag(triangles[i]->eval() * 1.E2);
    };

    xmin = 4.2;
    xmax = 4.4;
    ymin = -1.;
    ymax = 5.;
    xlabel   = "#sqrt{#it{s}}  [GeV]";
    filename = "triangle_s.pdf";

    i = 0;
    plotter->AddEntry(N, re, {xmin, xmax}, "Real", PRINT);
    i = 1;
    plotter->AddDashedEntry(N, re, {xmin, xmax}, PRINT);

    i = 0;
    plotter->AddEntry(N, im, {xmin, xmax}, "Imaginary", PRINT);
    i = 1;
    plotter->AddDashedEntry(N, im, {xmin, xmax}, PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.7,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    plotter->ClearData();
    external = {M_Y4260, M_ZC3900, M_PION};

    auto re2 = [&](double E)
    {
        external[1] = E;
        triangles[i]->set_external_masses(external);
        return real(triangles[i]->eval() * 1.E2);
    };
    auto im2 = [&](double E)
    {
        external[1] = E;
        triangles[i]->set_external_masses(external);
        return imag(triangles[i]->eval() * 1.E2);
    };

    xmin = 3.8;
    xmax = 4.1;
    ymin = 0.;
    ymax = 1.5;
    xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";
    filename = "triangle_sab.pdf";

    i = 0;
    plotter->AddEntry(N, re2, {xmin, xmax}, "Real", PRINT);
    i = 1;
    plotter->AddDashedEntry(N, re2, {xmin, xmax}, PRINT);

    i = 0;
    plotter->AddEntry(N, im2, {xmin, xmax}, "Imaginary", PRINT);
    i = 1;
    plotter->AddDashedEntry(N, im2, {xmin, xmax}, PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.7,0.7);
    plotter->SetLegendOffset(0.4, 0.07);
    plotter->Plot(filename);

    delete plotter;
};  