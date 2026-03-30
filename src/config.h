/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file defines how to use the config file
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CONFIG_H
#define SMASH_ROOT_ANALYSIS_CONFIG_H

#include <string>
#include <unordered_set>

#include "./constants.h"



struct Config {
  // Verbosity
  bool verbose = true;

  
  
  // Directories:
  int start_directory = 0;
  int number_of_directories = 10;


  
  // Analyses:
  // Flow
  bool flow_basic = false;


  
  // Multiplicity
  bool multiplicity = false; 
  // multiplicity analysis parameters
  std::unordered_set<int> multiplicity_excluded_species{};
  bool multiplicity_FXT_frame = false;
  bool multiplicity_default_cuts = false;
  double multiplicity_eta_min{};
  double multiplicity_eta_max{};
  double multiplicity_pT_min{};
  double multiplicity_pT_max{};
  std::vector<double> centrality_class_edges =
    {0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};

  void default_multiplicity_cuts() {
    if (multiplicity_FXT_frame) {
      multiplicity_eta_min = Eta_Multiplicity_FXT_Min;
      multiplicity_eta_max = Eta_Multiplicity_FXT_Max;
    } else {
      multiplicity_eta_min = Eta_Multiplicity_Min;
      multiplicity_eta_max = Eta_Multiplicity_Max;
    }

    multiplicity_pT_min = PT_Multiplicity_Min;
    multiplicity_pT_max = PT_Multiplicity_Max;
  }
  
  

  // Yields
  bool yields = false;
  // Assume no low pT cut
  double yields_proton_pT_min = PT_Baryon_Yield_Min;
  double yields_lambda_pT_min = PT_Baryon_Yield_Min;
  double yields_pi_pT_min = PT_Meson_Yield_Min;
  double yields_kaon_pT_min = PT_Meson_Yield_Min;
  double yields_phi_pT_min = PT_Meson_Yield_Min;

  

  // Read from file
  void load(const std::string& filename);
};



#endif
