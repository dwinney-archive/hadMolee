// Wrapper class to call box loop function in LoopTools
// This needs to be compiled as its own seperate library because it needs 
// to be linked with gFortran
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "LT_box.hpp"

complex<double> LT_box::eval()
{
    complex<double> boxx;

    ltini();
    if ( (_wA + _wB + _wC + _wD) < EPS ) {boxx = D0 (_emA2, _emD2, _emC2, _emB2, _t, _s, _imB2,             _imA2,             _imD2,             _imC2);}
    else                                 {boxx = D0C(_emA2, _emD2, _emC2, _emB2, _t, _s, _imB2-XI*_imB*_wB, _imA2-XI*_imA*_wA, _imD2-XI*_imD*_wD, _imC2-XI*_imC*_wC);}
    ltexi();
    
    return - boxx / pow(4.*PI, 2.);
};