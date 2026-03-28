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

#include "./constants.h"



struct Config {
  // Directories:
  int start_directory = 0;
  int number_of_directories = 10;

  // Analyses:
  // Multiplicity
  bool multiplicity = false;

  // Yields
  bool yields = false;
  // Assume no low pT cut
  double yields_proton_pT_min = PT_Baryon_Yield_Min;
  double yields_lambda_pT_min = PT_Baryon_Yield_Min;
  double yields_pi_pT_min = PT_Meson_Yield_Min;
  double yields_kaon_pT_min = PT_Meson_Yield_Min;
  double yields_phi_pT_min = PT_Meson_Yield_Min;

  // Flow
  bool flow_basic = false;

  // Read from file
  void load(const std::string& filename);
};



#endif
