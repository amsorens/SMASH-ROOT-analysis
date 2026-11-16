/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains helper structs and functions used in various parts of the code
// (config, multiple analyses, etc.).
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_SHARED_HELPERS_H
#define SMASH_ROOT_ANALYSIS_SHARED_HELPERS_H

#include "./basic_functions.h"
#include "./fourvector.h"
#include "./range.h"



/////////////////////////////////////////////
// Multiplicity

struct MultiplicityCuts {
  // list of particle pdg codes to exclude from multiplicity
  std::unordered_set<int> pdgs_to_exclude;
  // whether the multiplicity is computed in the FXT frame
  bool FXT_frame;
  // considered pseudorapidity range
  Range<double> eta;
  // considered transverse momentum pT range
  Range<double> pT;
};

bool passes_multiplicity_cuts(const MultiplicityCuts& multiplicity_cuts, int charge,
			      int pdgcode, FourVector p_mu, double v_beam);



#endif
