#include "triangle.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

using namespace std;

void triangles_plot()
{   
    array<double,3> internal = {M_DSTAR, M_D1,  M_D};
    array<double,3> external = {M_Y,     M_Z,   M_PION};

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

    // double xmin = 4.2;
    // double xmax = 4.4;
    // double ymin = -1.;
    // double ymax = 5.;
    // string xlabel   = "#sqrt{#it{s}}  [GeV]";
    
    double xmin = 3.8;
    double xmax = 4.1;
    double ymin = 0.;
    double ymax = 1.5;
    string xlabel   = "#sqrt{#sigma_{#it{D}*#it{D}}}  [GeV]";

    string ylabel   = "#it{T}_{#it{D}*#it{D}_{1}#it{D}} #times 100";
    string filename = "tri_test.pdf";

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    int i;
    auto re = [&](double E)
    {
        external[1] = E;
        triangles[i]->set_external_masses(external);
        return real(triangles[i]->eval() * 1.E2);
    };
    auto im = [&](double E)
    {
        external[1] = E;
        triangles[i]->set_external_masses(external);
        return imag(triangles[i]->eval() * 1.E2);
    };

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
    plotter->RemoveLogo();
    plotter->Plot(filename);

    delete plotter;
};  