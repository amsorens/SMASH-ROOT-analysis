/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains basic constants used;
// all constants from the constants file use capital letters.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CONSTANTS_H
#define SMASH_ROOT_ANALYSIS_CONSTANTS_H



//////////////////////////////////////////////////////////////////////////////////////////
// File management

// This analysis suite assumes the following structure: a directory "SMASH_results" (can
// be any other name) contains a directory called "data" (note: **must** be that name) and
// "analysis" (can be any other name). The "data" directory must contain results of SMASH
// simulations in directories named with integer numbers, which **must** start from zero
// (note that gaps in numbering are allowed, so that {"0", "1", "2", "3", "4", ...} and
// {"0", "3", "4", "5", "27", ...} are both fine). The analysis suite is copied into
// "analysis" as "SMASH_ROOT_analysis" (can be any other name). In "SMASH_ROOT_analysis",
// there is a "src" directory with the source code, and the user creates a "build"
// directory (can be any other name) where the analysis suite is compiled. Thus we have
// SMASH_results/                   # Top-level directory (can be any name)
// ├── data/                        # Must be named exactly "data"
// │   ├── 0/                       # Integer-named subdirectories (must start from 0)
// │   ├── 3/
// │   ├── 4/
// │   └── ...
// └── analysis/                    # Can be any name
//     └── SMASH_ROOT_analysis/     # Analysis suite (can be any name)
//         ├── src/                 # Source code directory (fixed name)
//         └── build/               # User-created build directory (can be any name)
// The target directory for the analysis suite output is the "analysis" directory, which
// is two levels above the "build" directory.
const char Target_Directory[256] = "../../";

// The size of char arrays used in file names etc.
const size_t Char_Array_Size = 256;

///////////////////////////////////////////
// Maximum number of particles in a single ROOT entry, needed to populate arrays of
// particle momenta, positions, etc. Generally, the size of these arrays should be equal
// to or larger than the number of particles. In practice, we don't know that number ahead
// of time. The default SMASH ROOT output has a limit of 500,000 particles per ROOT entry,
// and it works well for most applications. Therefore, we set an extremely conservative
// limit of 5,000,000. For other applications, where the number of particles per entry is
// even larger, the code will crash with an error message advising the user to increase
// Max_Particles.
const int Max_Particles = 5000000;



//////////////////////////////////////////////////////////////////////////////////////////
// Nuclear properties

// mN needs to be 0.938 as things are parameterized in SMASH with that value (?!)
const double Nucleon_Mass = 0.938; // in GeV



//////////////////////////////////////////////////////////////////////////////////////////
// Experimental cuts

/////////////////////////////////////////////
// STAR

// STAR BES I & II collider pseudorapidity and pT cuts for multiplicity
const double Eta_Multiplicity_Min = -1.0;
const double Eta_Multiplicity_Max = 1.0;
const double PT_Multiplicity_Min = 0.2; // [GeV]
const double PT_Multiplicity_Max = 50.0; // [GeV]
// STAR FXT pseudorapidity cuts for multiplicity (pT cuts the same as collider)
const double Eta_Multiplicity_FXT_Min = 0.05;
const double Eta_Multiplicity_FXT_Max = 2.0;

// STAR BES I & II collider pT cuts for proton cumulants
const double PT_Cumulants_Min = 0.4; // [GeV]
const double PT_Cumulants_Max = 2.0; // [GeV]

// STAR BES I & II collider pT cuts for yields and freeze-out parameters (CHECK!!!)
const double PT_Meson_Yield_Min = 0.2; // [GeV]
const double PT_Baryon_Yield_Min = 0.4; // [GeV]

// STAR BES FXT pT cuts for flow
const double STAR_FXT_proton_pT_min = 0.4; // [GeV]
const double STAR_FXT_proton_pT_max = 2.0; // [GeV]
const double STAR_FXT_deuteron_pT_min = 0.8; // [GeV]
const double STAR_FXT_deuteron_pT_max = 2.0; // [GeV]
const double STAR_FXT_lambda_pT_min = 0.4; // [GeV]
const double STAR_FXT_lambda_pT_max = 2.0; // [GeV]
const double STAR_FXT_pion_pT_min = 0.2; // [GeV]
const double STAR_FXT_pion_pT_max = 1.6; // [GeV]
const double STAR_FXT_kaon_pT_min = 0.4; // [GeV]
const double STAR_FXT_kaon_pT_max = 1.6; // [GeV]

/////////////////////////////////////////////
// HADES

// HADES pT cuts for flow
const double HADES_proton_pT_min = 0.2; // [GeV]
const double HADES_proton_pT_max = 2.0; // [GeV]
const double HADES_deuteron_pT_min = 0.8; // [GeV]
const double HADES_deuteron_pT_max = 2.0; // [GeV]
const double HADES_lambda_pT_min = 0.4; // [GeV]
const double HADES_lambda_pT_max = 2.0; // [GeV]
const double HADES_pion_pT_min = 0.2; // [GeV]
const double HADES_pion_pT_max = 1.6; // [GeV]
const double HADES_kaon_pT_min = 0.4; // [GeV]
const double HADES_kaon_pT_max = 1.6; // [GeV]

/////////////////////////////////////////////
// FOPI

// FOPI pT cuts for flow
// all of these are uncertain -- DOUBLE CHECK for a serious study
const double FOPI_proton_pT_min = 0.375; // [GeV]
const double FOPI_proton_pT_max = 10.0; // [GeV], they don't really use a max
const double FOPI_deuteron_pT_min = 0.375; // [GeV]
const double FOPI_deuteron_pT_max = 10.0; // [GeV], they don't really use a max
const double FOPI_lambda_pT_min = 0.375; // [GeV]
const double FOPI_lambda_pT_max = 10.0; // [GeV], they don't really use a max
const double FOPI_pion_pT_min = 0.2; // [GeV]
const double FOPI_pion_pT_max = 10.0; // [GeV], they don't really use a max
const double FOPI_kaon_pT_min = 0.2; // [GeV]
const double FOPI_kaon_pT_max = 10.0; // [GeV], they don't really use a max



#endif
