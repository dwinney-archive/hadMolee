// Wrapper class to call triangle loop function in LoopTools
// This needs to be compiled as its own seperate library because it needs 
// to be linked with gFortran
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "LT_triangle.hpp"

complex<double> LT_triangle::eval()
{
    complex<double> tri;

    ltini();
    if ( (_wA + _wB + _wC) > EPS ) tri = C0 (_emA2, _emB2, _emC2, _imB2, _imC2, _imA2);
    else                           tri = C0C(_emA2, _emB2, _emC2, _imB2 - XI*_imB*_wB, _imC2 - XI*_imC*_wC, _imA2 - XI*_imA*_wA);
    ltexi();
    
    return - tri / pow(4.*PI, 2.);
};