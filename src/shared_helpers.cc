/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains functions used in various parts of the code (config, multiple
// analyses, etc.).
//////////////////////////////////////////////////////////////////////////////////////////

#include "./shared_helpers.h"



/////////////////////////////////////////////
// Multiplicity

bool passes_multiplicity_cuts(const MultiplicityCuts& multiplicity_cuts, int charge,
			      int pdgcode, FourVector p_mu, double v_beam) {
  // Only charged particles
  if (charge == 0) {
    return false;
  }
  
  // Skip excluded particles
  if (multiplicity_cuts.pdgs_to_exclude.count(pdgcode)) {
    return false;
  }
		    
  // Apply the pT cut
  const double pT = transverse_momentum(p_mu.x1(), p_mu.x2());
  if (!(pT > multiplicity_cuts.pT.min() && pT < multiplicity_cuts.pT.max())) {
    return false;
  }
  
  // Apply the pseudorapidity cut
  if (multiplicity_cuts.FXT_frame) {
    // Boost to the fixed-target frame
    p_mu.boost_fourvector(0.0, 0.0, -v_beam);
  }      
  // Calculate pseudorapidity
  const double eta = pseudorapidity(p_mu.x1(), p_mu.x2(), p_mu.x3());	      
  return eta > multiplicity_cuts.eta.min() && eta < multiplicity_cuts.eta.max();
}



