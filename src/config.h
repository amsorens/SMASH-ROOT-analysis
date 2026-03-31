/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file defines how to use the config file
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CONFIG_H
#define SMASH_ROOT_ANALYSIS_CONFIG_H

#include <set>
#include <string>
#include <unordered_set>

#include "./constants.h"



struct Config {
  // Verbosity
  bool verbose = true;

  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Directories:
  int start_directory = 0;
  int number_of_directories = 10;


  ////////////////////////////////////////////////////////////////////////////////////////
  // Analyses:

  ///////////////////////////////////////////
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
  
  

  ///////////////////////////////////////////
  // Yields
  bool yields = false;
  // yields analysis parameters
  double yields_proton_pT_min = PT_Baryon_Yield_Min;
  double yields_lambda_pT_min = PT_Baryon_Yield_Min;
  double yields_pi_pT_min = PT_Meson_Yield_Min;
  double yields_kaon_pT_min = PT_Meson_Yield_Min;
  double yields_phi_pT_min = PT_Meson_Yield_Min;



  ///////////////////////////////////////////
  // Flow
  bool flow_basic = false;
  // Time evolution
  bool flow_only_at_final_output = true;  
  std::set<double> flow_output_times = {};
  // Rapidity coverage
  bool flow_scale_by_beam_rapidity = false;
  // Provide the number of rapidity bins as well as limits of the histogram range;
  // default values for 0.25 bin width, centered on -1.875, ... , -0.125, 0.125, ... :
  int flow_number_of_rapidity_bins = 16;
  double flow_y_min = -2.0;
  double flow_y_max = 2.0;
  // Coverage for quantities integrated at midrapidity:
  double flow_y_mid_min = -0.25;
  double flow_y_mid_max = 0.25;
  // Transverse momentum pT cuts
  bool flow_default_pT_cuts_HADES = false;
  bool flow_default_pT_cuts_STAR_FXT = false;
  double proton_pT_min{};
  double proton_pT_max{};
  double deuteron_pT_min{};
  double deuteron_pT_max{};
  double lambda_pT_min{};
  double lambda_pT_max{};
  double pion_pT_min{};
  double pion_pT_max{};
  double kaon_pT_min{};
  double kaon_pT_max{};

  void default_flow_pT_cuts() {
    if (flow_default_pT_cuts_HADES) {
      // Sanity check
      if (flow_default_pT_cuts_STAR_FXT) {
	throw std::runtime_error("Cannot ask for multiple default flow pT cuts.\n"
				 "Adjust the config file.");
      }      
      proton_pT_min = HADES_proton_pT_min;
      proton_pT_max = HADES_proton_pT_max;
      deuteron_pT_min = HADES_deuteron_pT_min;
      deuteron_pT_max = HADES_deuteron_pT_max;
      lambda_pT_min = HADES_lambda_pT_min;
      lambda_pT_max = HADES_lambda_pT_max;
      pion_pT_min = HADES_pion_pT_min;
      pion_pT_max = HADES_pion_pT_max;
      kaon_pT_min = HADES_kaon_pT_min;
      kaon_pT_max = HADES_kaon_pT_max;
    }

    if (flow_default_pT_cuts_STAR_FXT) {
      // Sanity check
      if (flow_default_pT_cuts_HADES) {
	throw std::runtime_error("Cannot ask for multiple default flow pT cuts.\n"
				 "Adjust the config file.");
      }
      proton_pT_min = STAR_FXT_proton_pT_min;
      proton_pT_max = STAR_FXT_proton_pT_max;
      deuteron_pT_min = STAR_FXT_deuteron_pT_min;
      deuteron_pT_max = STAR_FXT_deuteron_pT_max;
      lambda_pT_min = STAR_FXT_lambda_pT_min;
      lambda_pT_max = STAR_FXT_lambda_pT_max;
      pion_pT_min = STAR_FXT_pion_pT_min;
      pion_pT_max = STAR_FXT_pion_pT_max;
      kaon_pT_min = STAR_FXT_kaon_pT_min;
      kaon_pT_max = STAR_FXT_kaon_pT_max;
    }
  }


  
  // Read from file
  void load(const std::string& filename);
};



#endif
