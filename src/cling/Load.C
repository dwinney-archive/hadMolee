// Load script that links all the required libraries.
// Adapted from the installation of elSpectro
// [https://github.com/dglazier/elSpectro]
// by Derek Glazier and others
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// -----------------------------------------------------------------------------

void Load()
{
    TString LIB_EXT = gSystem->GetSoExt();
    TString BASE    = gSystem->Getenv("HADMOLEE");

    //----------------------------------------------------------------------
    // Core physics library
    TString PHYS_INC  = BASE + "/physics";

    TString CORE_LIB  = BASE + "/lib/libHADMOLEE." + LIB_EXT;
    TString CORE_INC  = BASE + "/src/core";

    if (!gSystem->AccessPathName(CORE_LIB.Data()))
    {
        gInterpreter->AddIncludePath( CORE_INC.Data());
        gInterpreter->AddIncludePath( PHYS_INC.Data());
        Int_t lib = gSystem->Load( CORE_LIB.Data());
    }
    else
    {
        Warning("Load", "hadMolee library not found! Path given: %s", CORE_LIB.Data());
    }

    //----------------------------------------------------------------------
    // Linked LT physics library

    // TString LT_INC  = BASE + "/src/LoopTools";
    // TString LT_LIB  = BASE + "/lib/libLOOPEXT." + LIB_EXT;

    // if (!gSystem->AccessPathName(LT_LIB.Data()))
    // {
    //     gInterpreter->AddIncludePath( LT_INC.Data());
    //     Int_t lib = gSystem->Load( LT_LIB.Data());
    // }
    // else
    // {
    //     Warning("Load", "LoopTools extention library not found! Path given: %s", LT_LIB.Data());
    // }
}