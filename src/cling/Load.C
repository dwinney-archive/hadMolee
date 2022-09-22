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
    TString BASE    = gSystem->Getenv("EE3BODY");

    //----------------------------------------------------------------------
    // Core physics library

    TString CORE_INC  = BASE;
            CORE_INC += "/include/core";
    TString CORE_LIB  = BASE;
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
    // Linked LT physics library

    TString LT_INC  = BASE;
            LT_INC += "/include/LoopTools";
    TString LT_LIB  = BASE;
            LT_LIB += "/lib/libLOOPEXT.";
            LT_LIB += LIB_EXT;

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