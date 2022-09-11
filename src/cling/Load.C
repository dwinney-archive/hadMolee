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

    //----------------------------------------------------------------------
    // Core physics library

    TString CORE_DIR  = gSystem->Getenv("EE3BODY");
    TString CORE_INC  = CORE_DIR;
            CORE_INC += "/include/";
    TString CORE_LIB  = CORE_DIR;
            CORE_LIB += "/lib/libEE3BODY.";
            CORE_LIB += LIB_EXT;

    if (!gSystem->AccessPathName(CORE_LIB.Data()))
    {
        gInterpreter->AddIncludePath( CORE_INC.Data());
        Int_t lib = gSystem->Load( CORE_LIB.Data());
    }
    else
    {
        Warning("Load", "ee3Body library not found! Path given: %s", CORE_LIB.Data());
    }

     //----------------------------------------------------------------------
    // Plotting library

    TString JPACSTYLE_DIR  = gSystem->Getenv("JPACSTYLE");
    TString JPACSTYLE_INC  = JPACSTYLE_DIR;
            JPACSTYLE_INC += "/include/";
    TString JPACSTYLE_LIB  = JPACSTYLE_DIR;
            JPACSTYLE_LIB += "/lib/libjpacStyle.";
            JPACSTYLE_LIB += LIB_EXT;

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