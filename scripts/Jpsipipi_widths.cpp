#include "kinematics.hpp"
#include "Y(4260).hpp"
#include "Jpsipipi_boxes.hpp"
#include "plotter.hpp"

using namespace hadMolee;

void Jpsipipi_widths()
{
    // Masses
    double ma = M_JPSI;
    double mb = M_PION;
    double mc = M_PION;

    kinematics kJpsipipi  = make_kinematics({ma, mb, mc});
    lineshape  Y          = make_lineshape<D1D_molecule>();

    amplitude  box_I      = make_amplitude<Jpsipipi_BoxI>(kJpsipipi, Y, "Box I");
};