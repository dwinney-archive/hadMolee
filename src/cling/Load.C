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

    TString CORE_LIB  = BASE + "/lib/libHADMOLEE." + LIB_EXT;
    TString CORE_INC  = BASE + "/include/core";
    TString PHYS_INC  = BASE + "/include/physics";

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

    TString LT_INC  = BASE + "/include/LoopTools";
    TString LT_LIB  = BASE + "/lib/libLOOPEXT." + LIB_EXT;

    if (!gSystem->AccessPathName(LT_LIB.Data()))
    {
        gInterpreter->AddIncludePath( LT_INC.Data());
        Int_t lib = gSystem->Load( LT_LIB.Data());
    }
    else
    {
        Warning("Load", "LoopTools extention library not found! Path given: %s", LT_LIB.Data());
    }

    //----------------------------------------------------------------------
    // Need to add cubature to apth
    TString CUBE_INC  = gSystem->Getenv("CUBATURE");
    gInterpreter->AddIncludePath( CUBE_INC.Data());

     //----------------------------------------------------------------------
    // Plotting library

    TString JPACSTYLE_DIR  = gSystem->Getenv("JPACSTYLE");
    TString JPACSTYLE_INC  = JPACSTYLE_DIR + "/include/";
    TString JPACSTYLE_LIB  = JPACSTYLE_DIR + "/lib/libjpacStyle." + LIB_EXT;

    if (!gSystem->AccessPathName(JPACSTYLE_LIB.Data()))
    {
        gInterpreter->AddIncludePath( JPACSTYLE_INC.Data());
        Int_t stylib = gSystem->Load( JPACSTYLE_LIB.Data());
    }
    else
    {
        Warning("Load", "jpacStyle library not found! Path given: %s", JPACSTYLE_LIB.Data());
    }
}