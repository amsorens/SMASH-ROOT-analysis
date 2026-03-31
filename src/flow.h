/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate collective flow.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_FLOW_H
#define SMASH_ROOT_ANALYSIS_FLOW_H

#include <iostream>
#include <memory>
#include <set>

#include <TROOT.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>

#include "./basic_functions.h"
#include "./config.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"



inline void initialize_tgraph
  (std::unique_ptr<TGraphErrors>& tgraph_to_initialize,
   const char* basic_name, const char* prefix) {
  char graph_name[Char_Array_Size];
  snprintf(graph_name, Char_Array_Size, "%s%s", prefix, basic_name);

  tgraph_to_initialize = std::make_unique<TGraphErrors>();
  tgraph_to_initialize->SetName(graph_name);
  tgraph_to_initialize->SetTitle(graph_name);
}



inline void populate_tgraphs
  (std::vector<std::unique_ptr<TGraphErrors>>& vector_of_tgraphs_to_populate,
   const char* basic_name, const char* prefix) {
  char graph_name[Char_Array_Size];
  snprintf(graph_name, Char_Array_Size, "%s%s", prefix, basic_name);

  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(graph_name);
  graph->SetTitle(graph_name);

  vector_of_tgraphs_to_populate.push_back(std::move(graph));
}



class Flow {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Default constructor, "implemented" here in the .h file
  // We provide default values of the number of bins and spread in rapidity y for the
  // TProfile histograms, as well as the default value of bool for displaying stats
  Flow(double sqrts,
       //////////////////////////////////////
       // Time evolution of flow:
       bool only_flow_at_final_output,
       int number_of_time_steps,
       //////////////////////////////////////
       // Rapidity coverage:
       bool scale_by_beam_rapidity,
       int y_number_of_bins,	       
       double y_min,
       double y_max,
       //////////////////////////////////////
       // Coverage for quantities integrated at midrapidity:
       double y_mid_min,
       double y_mid_max,
       //////////////////////////////////////
       // Transverse momentum pT cuts:
       double proton_pT_min,
       double proton_pT_max,
       double deuteron_pT_min,
       double deuteron_pT_max,
       double lambda_pT_min,
       double lambda_pT_max,
       double pion_pT_min,
       double pion_pT_max,
       double kaon_pT_min,
       double kaon_pT_max,
       //////////////////////////////////////
       // Plotting parameters:
       bool show_stats = false,
       bool user_axes_ranges = true,
       // Plotting parameters; weird default values are updated in the constructor:
       double x_axis_lower = -100.0,
       double x_axis_upper  = -100.0,
       double v1_y_axis_lower = -100.0,
       double v1_y_axis_upper = - 100.0,
       double v2_y_axis_lower = -100.0,
       double v2_y_axis_upper = -100.0,
       double v3_y_axis_lower = -100.0,
       double v3_y_axis_upper = -100.0,
       double fitting_range = -100.0)
    // Explicit initialization of const terms (terms are initialized in the order of
    // declaration in the class definition, regardless of the order of appearance here):
    : sqrts_(sqrts),
      Ekin_( Ekin_from_sqrts(sqrts_) ),
      y_beam_( y_beam_cm(sqrts_) ),
      ///////////////////////////////////////
      // Time evolution of flow:
      only_flow_at_final_output_(only_flow_at_final_output),
      number_of_time_steps_(number_of_time_steps),
      ///////////////////////////////////////
      // Rapidity coverage:
      scale_by_beam_rapidity_(scale_by_beam_rapidity),
      y_number_of_bins_(y_number_of_bins),
      // note: using different values may affect which bins go into evaluation of
      // observables at midrapidity (e.g., |y|<0.5) etc.
      // TO DO: consider what to do if you scale by beam rapidity
      y_min_(y_min),
      y_max_(y_max),
      y_mid_min_(y_mid_min),
      y_mid_max_(y_mid_max),
      ///////////////////////////////////////
      // pT cuts:      
      proton_pT_min_(proton_pT_min),
      proton_pT_max_(proton_pT_max),
      deuteron_pT_min_(deuteron_pT_min),
      deuteron_pT_max_(deuteron_pT_max),
      lambda_pT_min_(lambda_pT_min),
      lambda_pT_max_(lambda_pT_max),
      pion_pT_min_(pion_pT_min),
      pion_pT_max_(pion_pT_max),
      kaon_pT_min_(kaon_pT_min),
      kaon_pT_max_(kaon_pT_max),
      ///////////////////////////////////////
      // Plotting parameters:
      user_axes_ranges_(user_axes_ranges),
      // use a lambda function: if axes ranges are not given, specify using y_beam_
      x_axis_lower_( ([&]() -> double {
	if (x_axis_lower == -100.0) { return -1.5 * y_beam_; }
	else { return x_axis_lower; }
      }) () ),
      x_axis_upper_( ([&]() -> double {
	if (x_axis_upper == -100.0) { return 1.5 * y_beam_; }
	else { return x_axis_upper; }
      }) () ),
      v1_y_axis_lower_( ([&]() -> double {
	if (v1_y_axis_lower == -100.0) { return -0.75; }
	else { return v1_y_axis_lower; }
      }) () ),
      v1_y_axis_upper_( ([&]() -> double {
	if (v1_y_axis_upper == -100.0) { return 0.75; }
	else { return v1_y_axis_upper; }
      }) () ),
      v2_y_axis_lower_( ([&]() -> double {
	if (v2_y_axis_lower == -100.0) { return -0.15; }
	else { return v2_y_axis_lower; }
      }) () ),
      v2_y_axis_upper_( ([&]() -> double {
	if (v2_y_axis_upper == -100.0) { return 0.15; }
	else { return v2_y_axis_upper; }
      }) () ),
      v3_y_axis_lower_( ([&]() -> double {
	if (v3_y_axis_lower == -100.0) { return -0.05; }
	else { return v3_y_axis_lower; }
      }) () ),
      v3_y_axis_upper_( ([&]() -> double {
	if (v3_y_axis_upper == -100.0) { return 0.05; }
	else { return v3_y_axis_upper; }
      }) () ),
      fitting_range_( ([&]() -> double {
	if (fitting_range == -100.0) { return 0.5 * y_beam_; }
	else { return fitting_range; }
      }) () ) {
    /////////////////////////////////////////
    // Sanity check for time evolution
    if ( !only_flow_at_final_output && (number_of_time_steps == 1) ) {
      throw std::runtime_error("You must provide number_of_time_steps > 1 to "
			       "calculate flow at all time steps.");
    }
    //////////////////////////////////////////////////////////////////////////////////////
    // Initialize all elements of TProfile vectors: one TProfile for each time step for a
    // given observable;
    for (int i = 0; i < number_of_time_steps_; i++) {
      ////////////////////////////////////////////////////////////////////////////////////
      // Initialize the yield histograms
      char h_N_proton_in_4pi_name[Char_Array_Size];
      snprintf(h_N_proton_in_4pi_name, Char_Array_Size, "h_N_proton_in_4pi_%d", i);
      h_N_proton_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_proton_in_4pi_name, h_N_proton_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_proton_in_4pi_[i]->Sumw2();   
      h_N_proton_in_4pi_[i]->SetStats(show_stats);

      char h_N_deuteron_in_4pi_name[Char_Array_Size];
      snprintf(h_N_deuteron_in_4pi_name, Char_Array_Size, "h_N_deuteron_in_4pi_%d", i);
      h_N_deuteron_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_deuteron_in_4pi_name, h_N_deuteron_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_deuteron_in_4pi_[i]->Sumw2();   
      h_N_deuteron_in_4pi_[i]->SetStats(show_stats);

      char h_N_lambda_in_4pi_name[Char_Array_Size];
      snprintf(h_N_lambda_in_4pi_name, Char_Array_Size, "h_N_lambda_in_4pi_%d", i);
      h_N_lambda_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_lambda_in_4pi_name, h_N_lambda_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_lambda_in_4pi_[i]->Sumw2();   
      h_N_lambda_in_4pi_[i]->SetStats(show_stats);

      char h_N_pi_plus_in_4pi_name[Char_Array_Size];
      snprintf(h_N_pi_plus_in_4pi_name, Char_Array_Size, "h_N_pi_plus_in_4pi_%d", i);
      h_N_pi_plus_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_pi_plus_in_4pi_name, h_N_pi_plus_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_pi_plus_in_4pi_[i]->Sumw2();   
      h_N_pi_plus_in_4pi_[i]->SetStats(show_stats);

      char h_N_pi_minus_in_4pi_name[Char_Array_Size];
      snprintf(h_N_pi_minus_in_4pi_name, Char_Array_Size, "h_N_pi_minus_in_4pi_%d", i);
      h_N_pi_minus_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_pi_minus_in_4pi_name, h_N_pi_minus_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_pi_minus_in_4pi_[i]->Sumw2();   
      h_N_pi_minus_in_4pi_[i]->SetStats(show_stats);

      char h_N_kaon_plus_in_4pi_name[Char_Array_Size];
      snprintf(h_N_kaon_plus_in_4pi_name, Char_Array_Size, "h_N_kaon_plus_in_4pi_%d", i);
      h_N_kaon_plus_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_kaon_plus_in_4pi_name, h_N_kaon_plus_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_kaon_plus_in_4pi_[i]->Sumw2();   
      h_N_kaon_plus_in_4pi_[i]->SetStats(show_stats);

      char h_N_kaon_minus_in_4pi_name[Char_Array_Size];
      snprintf(h_N_kaon_minus_in_4pi_name, Char_Array_Size, "h_N_kaon_minus_in_4pi_%d", i);
      h_N_kaon_minus_in_4pi_.emplace_back
	(std::make_unique<TH1D>(h_N_kaon_minus_in_4pi_name, h_N_kaon_minus_in_4pi_name,
				y_number_of_bins_, y_min_, y_max_));
      h_N_kaon_minus_in_4pi_[i]->Sumw2();   
      h_N_kaon_minus_in_4pi_[i]->SetStats(show_stats);

      // Initialize the yield histograms at midrapidity; the number of bins is hardcoded:
      char h_N_proton_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_proton_at_mid_y_name, Char_Array_Size, "h_N_proton_at_mid_y_%d", i);
      h_N_proton_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_proton_at_mid_y_name, h_N_proton_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_proton_at_mid_y_[i]->Sumw2();   
      h_N_proton_at_mid_y_[i]->SetStats(show_stats);

      char h_N_deuteron_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_deuteron_at_mid_y_name, Char_Array_Size, "h_N_deuteron_at_mid_y_%d", i);
      h_N_deuteron_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_deuteron_at_mid_y_name, h_N_deuteron_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_deuteron_at_mid_y_[i]->Sumw2();   
      h_N_deuteron_at_mid_y_[i]->SetStats(show_stats);

      char h_N_lambda_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_lambda_at_mid_y_name, Char_Array_Size, "h_N_lambda_at_mid_y_%d", i);
      h_N_lambda_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_lambda_at_mid_y_name, h_N_lambda_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_lambda_at_mid_y_[i]->Sumw2();   
      h_N_lambda_at_mid_y_[i]->SetStats(show_stats);

      char h_N_pi_plus_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_pi_plus_at_mid_y_name, Char_Array_Size, "h_N_pi_plus_at_mid_y_%d", i);
      h_N_pi_plus_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_pi_plus_at_mid_y_name, h_N_pi_plus_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_pi_plus_at_mid_y_[i]->Sumw2();   
      h_N_pi_plus_at_mid_y_[i]->SetStats(show_stats);

      char h_N_pi_minus_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_pi_minus_at_mid_y_name, Char_Array_Size, "h_N_pi_minus_at_mid_y_%d", i);
      h_N_pi_minus_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_pi_minus_at_mid_y_name, h_N_pi_minus_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_pi_minus_at_mid_y_[i]->Sumw2();   
      h_N_pi_minus_at_mid_y_[i]->SetStats(show_stats);

      char h_N_kaon_plus_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_kaon_plus_at_mid_y_name, Char_Array_Size, "h_N_kaon_plus_at_mid_y_%d", i);
      h_N_kaon_plus_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_kaon_plus_at_mid_y_name, h_N_kaon_plus_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_kaon_plus_at_mid_y_[i]->Sumw2();   
      h_N_kaon_plus_at_mid_y_[i]->SetStats(show_stats);

      char h_N_kaon_minus_at_mid_y_name[Char_Array_Size];
      snprintf(h_N_kaon_minus_at_mid_y_name, Char_Array_Size, "h_N_kaon_minus_at_mid_y_%d", i);
      h_N_kaon_minus_at_mid_y_.emplace_back
	(std::make_unique<TH1D>(h_N_kaon_minus_at_mid_y_name, h_N_kaon_minus_at_mid_y_name,
				10, y_mid_min_, y_mid_max_));
      h_N_kaon_minus_at_mid_y_[i]->Sumw2();   
      h_N_kaon_minus_at_mid_y_[i]->SetStats(show_stats);
      
      ////////////////////////////////////////////////////////////////////////////////////
      // Initialize the v1 TProfiles
      char p_v1_proton_basic_name[Char_Array_Size];
      snprintf(p_v1_proton_basic_name, Char_Array_Size, "p_v1_proton_basic_%d", i);
      p_v1_proton_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_proton_basic_name, p_v1_proton_basic_name,
				    // number_of_bins of bins from y_min to y_max, filled
				    // with values from -1 to 1
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_proton_basic_[i]->Sumw2();   
      p_v1_proton_basic_[i]->SetStats(show_stats);

      char p_v1_deuteron_basic_name[Char_Array_Size];
      snprintf(p_v1_deuteron_basic_name, Char_Array_Size, "p_v1_deuteron_basic_%d", i);
      p_v1_deuteron_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_deuteron_basic_name, p_v1_deuteron_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_deuteron_basic_[i]->Sumw2();   
      p_v1_deuteron_basic_[i]->SetStats(show_stats);

      char p_v1_lambda_basic_name[Char_Array_Size];
      snprintf(p_v1_lambda_basic_name, Char_Array_Size, "p_v1_lambda_basic_%d", i);
      p_v1_lambda_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_lambda_basic_name, p_v1_lambda_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_lambda_basic_[i]->Sumw2();   
      p_v1_lambda_basic_[i]->SetStats(show_stats);

      char p_v1_pi_plus_basic_name[Char_Array_Size];
      snprintf(p_v1_pi_plus_basic_name, Char_Array_Size, "p_v1_pi_plus_basic_%d", i);
      p_v1_pi_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_pi_plus_basic_name, p_v1_pi_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_pi_plus_basic_[i]->Sumw2();
      p_v1_pi_plus_basic_[i]->SetStats(show_stats);

      char p_v1_pi_minus_basic_name[Char_Array_Size];
      snprintf(p_v1_pi_minus_basic_name, Char_Array_Size, "p_v1_pi_minus_basic_%d", i);
      p_v1_pi_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_pi_minus_basic_name, p_v1_pi_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_pi_minus_basic_[i]->Sumw2();
      p_v1_pi_minus_basic_[i]->SetStats(show_stats);

      char p_v1_kaon_plus_basic_name[Char_Array_Size];
      snprintf(p_v1_kaon_plus_basic_name, Char_Array_Size, "p_v1_kaon_plus_basic_%d", i);
      p_v1_kaon_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_kaon_plus_basic_name, p_v1_kaon_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_kaon_plus_basic_[i]->Sumw2();
      p_v1_kaon_plus_basic_[i]->SetStats(show_stats);

      char p_v1_kaon_minus_basic_name[Char_Array_Size];
      snprintf(p_v1_kaon_minus_basic_name, Char_Array_Size, "p_v1_kaon_minus_basic_%d", i);
      p_v1_kaon_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v1_kaon_minus_basic_name, p_v1_kaon_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v1_kaon_minus_basic_[i]->Sumw2();
      p_v1_kaon_minus_basic_[i]->SetStats(show_stats);

      ////////////////////////////////////////////////////////////////////////////////////
      // Initialize the v2 TProfiles
      char p_v2_proton_basic_name[Char_Array_Size];
      snprintf(p_v2_proton_basic_name, Char_Array_Size, "p_v2_proton_basic_%d", i);
      p_v2_proton_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_proton_basic_name, p_v2_proton_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_proton_basic_[i]->Sumw2();   
      p_v2_proton_basic_[i]->SetStats(show_stats);

      char p_v2_deuteron_basic_name[Char_Array_Size];
      snprintf(p_v2_deuteron_basic_name, Char_Array_Size, "p_v2_deuteron_basic_%d", i);
      p_v2_deuteron_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_deuteron_basic_name, p_v2_deuteron_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_deuteron_basic_[i]->Sumw2();   
      p_v2_deuteron_basic_[i]->SetStats(show_stats);

      char p_v2_lambda_basic_name[Char_Array_Size];
      snprintf(p_v2_lambda_basic_name, Char_Array_Size, "p_v2_lambda_basic_%d", i);
      p_v2_lambda_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_lambda_basic_name, p_v2_lambda_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_lambda_basic_[i]->Sumw2();   
      p_v2_lambda_basic_[i]->SetStats(show_stats);

      char p_v2_pi_plus_basic_name[Char_Array_Size];
      snprintf(p_v2_pi_plus_basic_name, Char_Array_Size, "p_v2_pi_plus_basic_%d", i);
      p_v2_pi_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_pi_plus_basic_name, p_v2_pi_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_pi_plus_basic_[i]->Sumw2();   
      p_v2_pi_plus_basic_[i]->SetStats(show_stats);

      char p_v2_pi_minus_basic_name[Char_Array_Size];
      snprintf(p_v2_pi_minus_basic_name, Char_Array_Size, "p_v2_pi_minus_basic_%d", i);
      p_v2_pi_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_pi_minus_basic_name, p_v2_pi_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_pi_minus_basic_[i]->Sumw2();   
      p_v2_pi_minus_basic_[i]->SetStats(show_stats);

      char p_v2_kaon_plus_basic_name[Char_Array_Size];
      snprintf(p_v2_kaon_plus_basic_name, Char_Array_Size, "p_v2_kaon_plus_basic_%d", i);
      p_v2_kaon_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_kaon_plus_basic_name, p_v2_kaon_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_kaon_plus_basic_[i]->Sumw2();   
      p_v2_kaon_plus_basic_[i]->SetStats(show_stats);

      char p_v2_kaon_minus_basic_name[Char_Array_Size];
      snprintf(p_v2_kaon_minus_basic_name, Char_Array_Size, "p_v2_kaon_minus_basic_%d", i);
      p_v2_kaon_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v2_kaon_minus_basic_name, p_v2_kaon_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v2_kaon_minus_basic_[i]->Sumw2();   
      p_v2_kaon_minus_basic_[i]->SetStats(show_stats);

      ////////////////////////////////////////////////////////////////////////////////////
      // Initialize the v3 TProfiles
      char p_v3_proton_basic_name[Char_Array_Size];
      snprintf(p_v3_proton_basic_name, Char_Array_Size, "p_v3_proton_basic_%d", i);
      p_v3_proton_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_proton_basic_name, p_v3_proton_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_proton_basic_[i]->Sumw2();   
      p_v3_proton_basic_[i]->SetStats(show_stats);

      char p_v3_deuteron_basic_name[Char_Array_Size];
      snprintf(p_v3_deuteron_basic_name, Char_Array_Size, "p_v3_deuteron_basic_%d", i);
      p_v3_deuteron_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_deuteron_basic_name, p_v3_deuteron_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_deuteron_basic_[i]->Sumw2();   
      p_v3_deuteron_basic_[i]->SetStats(show_stats);

      char p_v3_lambda_basic_name[Char_Array_Size];
      snprintf(p_v3_lambda_basic_name, Char_Array_Size, "p_v3_lambda_basic_%d", i);
      p_v3_lambda_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_lambda_basic_name, p_v3_lambda_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_lambda_basic_[i]->Sumw2();   
      p_v3_lambda_basic_[i]->SetStats(show_stats);

      char p_v3_pi_plus_basic_name[Char_Array_Size];
      snprintf(p_v3_pi_plus_basic_name, Char_Array_Size, "p_v3_pi_plus_basic_%d", i);
      p_v3_pi_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_pi_plus_basic_name, p_v3_pi_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_pi_plus_basic_[i]->Sumw2();   
      p_v3_pi_plus_basic_[i]->SetStats(show_stats);

      char p_v3_pi_minus_basic_name[Char_Array_Size];
      snprintf(p_v3_pi_minus_basic_name, Char_Array_Size, "p_v3_pi_minus_basic_%d", i);
      p_v3_pi_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_pi_minus_basic_name, p_v3_pi_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_pi_minus_basic_[i]->Sumw2();   
      p_v3_pi_minus_basic_[i]->SetStats(show_stats);

      char p_v3_kaon_plus_basic_name[Char_Array_Size];
      snprintf(p_v3_kaon_plus_basic_name, Char_Array_Size, "p_v3_kaon_plus_basic_%d", i);
      p_v3_kaon_plus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_kaon_plus_basic_name, p_v3_kaon_plus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_kaon_plus_basic_[i]->Sumw2();   
      p_v3_kaon_plus_basic_[i]->SetStats(show_stats);

      char p_v3_kaon_minus_basic_name[Char_Array_Size];
      snprintf(p_v3_kaon_minus_basic_name, Char_Array_Size, "p_v3_kaon_minus_basic_%d", i);
      p_v3_kaon_minus_basic_.emplace_back
	(std::make_unique<TProfile>(p_v3_kaon_minus_basic_name, p_v3_kaon_minus_basic_name,
				    y_number_of_bins_, y_min_, y_max_, -1.0, 1.0));
      p_v3_kaon_minus_basic_[i]->Sumw2();   
      p_v3_kaon_minus_basic_[i]->SetStats(show_stats);
      ////////////////////////////////////////////////////////////////////////////////////
      // Initialize the integrated v1, v2, v3 TProfiles
      char p_v1_proton_basic_integrated_in_4pi_name[Char_Array_Size];
      snprintf(p_v1_proton_basic_integrated_in_4pi_name, Char_Array_Size,
	       "p_v1_proton_basic_integrated_in_4pi_%d", i);
      p_v1_proton_basic_integrated_in_4pi_.emplace_back
	(std::make_unique<TProfile>(p_v1_proton_basic_integrated_in_4pi_name,
				    p_v1_proton_basic_integrated_in_4pi_name,
				    // 1 bin from y_min to y_max, filled
				    // with values from -1 to 1
				    1, y_min_, y_max_, -1.0, 1.0));
      p_v1_proton_basic_integrated_in_4pi_[i]->Sumw2();   
      p_v1_proton_basic_integrated_in_4pi_[i]->SetStats(show_stats);

      char p_v1_proton_basic_integrated_at_mid_y_name[Char_Array_Size];
      snprintf(p_v1_proton_basic_integrated_at_mid_y_name, Char_Array_Size,
	       "p_v1_proton_basic_integrated_at_mid_y_%d", i);
      p_v1_proton_basic_integrated_at_mid_y_.emplace_back
	(std::make_unique<TProfile>(p_v1_proton_basic_integrated_at_mid_y_name,
				    p_v1_proton_basic_integrated_at_mid_y_name,
				    // 1 bin from y_min to y_max, filled
				    // with values from -1 to 1
				    1, y_mid_min_, y_mid_max_, -1.0, 1.0));
      p_v1_proton_basic_integrated_at_mid_y_[i]->Sumw2();   
      p_v1_proton_basic_integrated_at_mid_y_[i]->SetStats(show_stats);

      char p_v2_proton_basic_integrated_in_4pi_name[Char_Array_Size];
      snprintf(p_v2_proton_basic_integrated_in_4pi_name, Char_Array_Size,
	       "p_v2_proton_basic_integrated_in_4pi_%d", i);
      p_v2_proton_basic_integrated_in_4pi_.emplace_back
	(std::make_unique<TProfile>(p_v2_proton_basic_integrated_in_4pi_name,
				    p_v2_proton_basic_integrated_in_4pi_name,
				    1, y_min_, y_max_, -1.0, 1.0));
      p_v2_proton_basic_integrated_in_4pi_[i]->Sumw2();   
      p_v2_proton_basic_integrated_in_4pi_[i]->SetStats(show_stats);

      char p_v2_proton_basic_integrated_at_mid_y_name[Char_Array_Size];
      snprintf(p_v2_proton_basic_integrated_at_mid_y_name, Char_Array_Size,
	       "p_v2_proton_basic_integrated_at_mid_y_%d", i);
      p_v2_proton_basic_integrated_at_mid_y_.emplace_back
	(std::make_unique<TProfile>(p_v2_proton_basic_integrated_at_mid_y_name,
				    p_v2_proton_basic_integrated_at_mid_y_name,
				    1, y_mid_min_, y_mid_max_, -1.0, 1.0));
      p_v2_proton_basic_integrated_at_mid_y_[i]->Sumw2();   
      p_v2_proton_basic_integrated_at_mid_y_[i]->SetStats(show_stats);

      char p_v3_proton_basic_integrated_in_4pi_name[Char_Array_Size];
      snprintf(p_v3_proton_basic_integrated_in_4pi_name, Char_Array_Size,
	       "p_v3_proton_basic_integrated_in_4pi_%d", i);
      p_v3_proton_basic_integrated_in_4pi_.emplace_back
	(std::make_unique<TProfile>(p_v3_proton_basic_integrated_in_4pi_name,
				    p_v3_proton_basic_integrated_in_4pi_name,
				    1, y_min_, y_max_, -1.0, 1.0));
      p_v3_proton_basic_integrated_in_4pi_[i]->Sumw2();   
      p_v3_proton_basic_integrated_in_4pi_[i]->SetStats(show_stats);

      char p_v3_proton_basic_integrated_at_mid_y_name[Char_Array_Size];
      snprintf(p_v3_proton_basic_integrated_at_mid_y_name, Char_Array_Size,
	       "p_v3_proton_basic_integrated_at_mid_y_%d", i);
      p_v3_proton_basic_integrated_at_mid_y_.emplace_back
	(std::make_unique<TProfile>(p_v3_proton_basic_integrated_at_mid_y_name,
				    p_v3_proton_basic_integrated_at_mid_y_name,
				    1, y_mid_min_, y_mid_max_, -1.0, 1.0));
      p_v3_proton_basic_integrated_at_mid_y_[i]->Sumw2();   
      p_v3_proton_basic_integrated_at_mid_y_[i]->SetStats(show_stats);
      
    }

    // TO DO: make vectors with all time steps also for other flow calculations

    //////////////////////////////////////////////////////////////////////////////////////
    // Initialize all TGraphErrors used for the time evolution of basic flow:
    // Binned in rapidity
    for (int i = 0; i < y_number_of_bins_; i++) {
      // Establish the bin center
      const double y_bin = y_min_ + (0.5 + i) * (y_max_ - y_min_)/y_number_of_bins_;
      //std::cout << "y_bin (center) = " << y_bin << std::endl;
      // Establish the basic name
      char basic_proton_name[Char_Array_Size];
      snprintf(basic_proton_name, Char_Array_Size, "proton_basic_y_%.3f", y_bin);
      // dvn_{coll}/dt:
      populate_tgraphs(g_dv1coll_dt_proton_basic_, basic_proton_name, "g_dv1colldt_");
      populate_tgraphs(g_dv2coll_dt_proton_basic_, basic_proton_name, "g_dv2colldt_");
      populate_tgraphs(g_dv3coll_dt_proton_basic_, basic_proton_name, "g_dv3colldt_");
      // dvn_{MF}/dt:
      populate_tgraphs(g_dv1MF_dt_proton_basic_, basic_proton_name, "g_dv1MFdt_");
      populate_tgraphs(g_dv2MF_dt_proton_basic_, basic_proton_name, "g_dv2MFdt_");
      populate_tgraphs(g_dv3MF_dt_proton_basic_, basic_proton_name, "g_dv3MFdt_");
      // dvntotal/dt:
      populate_tgraphs(g_dv1total_dt_proton_basic_, basic_proton_name, "g_dv1totaldt_");
      populate_tgraphs(g_dv2total_dt_proton_basic_, basic_proton_name, "g_dv2totaldt_");
      populate_tgraphs(g_dv3total_dt_proton_basic_, basic_proton_name, "g_dv3totaldt_");
      // vn_{coll}(t):
      populate_tgraphs(g_v1coll_proton_basic_, basic_proton_name, "g_v1coll_");
      populate_tgraphs(g_v2coll_proton_basic_, basic_proton_name, "g_v2coll_");
      populate_tgraphs(g_v3coll_proton_basic_, basic_proton_name, "g_v3coll_");;
      // vn_{MF}(t):
      populate_tgraphs(g_v1MF_proton_basic_, basic_proton_name, "g_MFcoll_");
      populate_tgraphs(g_v2MF_proton_basic_, basic_proton_name, "g_MFcoll_");
      populate_tgraphs(g_v3MF_proton_basic_, basic_proton_name, "g_MFcoll_");
      // vntotal(t):
      populate_tgraphs(g_v1total_proton_basic_, basic_proton_name, "g_totalcoll_");
      populate_tgraphs(g_v2total_proton_basic_, basic_proton_name, "g_totalcoll_");
      populate_tgraphs(g_v3total_proton_basic_, basic_proton_name, "g_totalcoll_");
    }
    // Integrated in rapidity
    // dvn_{coll}/dt:
    initialize_tgraph(g_dv1coll_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv1colldt_integrated_in_4pi_");
    initialize_tgraph(g_dv2coll_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv2colldt_integrated_in_4pi_");
    initialize_tgraph(g_dv3coll_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv3colldt_integrated_in_4pi_");
    initialize_tgraph(g_dv1coll_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv1colldt_integrated_at_mid_y_");
    initialize_tgraph(g_dv2coll_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv2colldt_integrated_at_mid_y_");
    initialize_tgraph(g_dv3coll_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv3colldt_integrated_at_mid_y_");
    // dvn_{MF}/dt:
    initialize_tgraph(g_dv1MF_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv1MFdt_integrated_in_4pi_");
    initialize_tgraph(g_dv2MF_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv2MFdt_integrated_in_4pi_");
    initialize_tgraph(g_dv3MF_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv3MFdt_integrated_in_4pi_");
    initialize_tgraph(g_dv1MF_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv1MFdt_integrated_at_mid_y_");
    initialize_tgraph(g_dv2MF_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv2MFdt_integrated_at_mid_y_");
    initialize_tgraph(g_dv3MF_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv3MFdt_integrated_at_mid_y_");
    // dvntotal/dt:
    initialize_tgraph(g_dv1total_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv1totaldt_integrated_in_4pi_");
    initialize_tgraph(g_dv2total_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv2totaldt_integrated_in_4pi_");
    initialize_tgraph(g_dv3total_dt_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_dv3totaldt_integrated_in_4pi_");
    initialize_tgraph(g_dv1total_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv1totaldt_integrated_at_mid_y_");
    initialize_tgraph(g_dv2total_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv2totaldt_integrated_at_mid_y_");
    initialize_tgraph(g_dv3total_dt_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_dv3totaldt_integrated_at_mid_y_");
    // vn_{coll}(t):
    initialize_tgraph(g_v1coll_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v1coll_integrated_in_4pi_");
    initialize_tgraph(g_v2coll_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v2coll_integrated_in_4pi_");
    initialize_tgraph(g_v3coll_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v3coll_integrated_in_4pi_");
    initialize_tgraph(g_v1coll_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v1coll_integrated_at_mid_y_");
    initialize_tgraph(g_v2coll_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v2coll_integrated_at_mid_y_");
    initialize_tgraph(g_v3coll_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v3coll_integrated_at_mid_y_");
    // vn_{MF}(t):
    initialize_tgraph(g_v1MF_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v1MF_integrated_in_4pi_");
    initialize_tgraph(g_v2MF_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v2MF_integrated_in_4pi_");
    initialize_tgraph(g_v3MF_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v3MF_integrated_in_4pi_");
    initialize_tgraph(g_v1MF_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v1MF_integrated_at_mid_y_");
    initialize_tgraph(g_v2MF_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v2MF_integrated_at_mid_y_");
    initialize_tgraph(g_v3MF_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v3MF_integrated_at_mid_y_");
    // vntotal(t):
    initialize_tgraph(g_v1total_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v1total_integrated_in_4pi_");
    initialize_tgraph(g_v2total_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v2total_integrated_in_4pi_");
    initialize_tgraph(g_v3total_proton_basic_integrated_in_4pi_,
		      "proton_basic_", "g_v3total_integrated_in_4pi_");
    initialize_tgraph(g_v1total_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v1total_integrated_at_mid_y_");
    initialize_tgraph(g_v2total_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v2total_integrated_at_mid_y_");
    initialize_tgraph(g_v3total_proton_basic_integrated_at_mid_y_,
		      "proton_basic_", "g_v3total_integrated_at_mid_y_");
  }

  // Default destructor
  virtual ~Flow() {}


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Member functions of the class  
  ////////////////////////////////////////////////////////////////////////////////////////

  // Function that fills dN/dy histograms with data
  void fill_N_vs_y
    (double p_0, double p_x, double p_y, double p_z,
     // pass by reference to avoid copy; const to avoid changes
     const std::unique_ptr<TH1D>& dN_dy_histogram,
     double low_pT_cut, double high_pT_cut);

  // Function that fills the v1 and v2 profile histograms with data
  void fill_v1_v2
    (double p_0, double p_x, double p_y, double p_z,
     // pass by reference to avoid copy; const to avoid changes
     const std::unique_ptr<TProfile>& v1_profile_histogram,
     const std::unique_ptr<TProfile>& v2_profile_histogram,
     double low_pT_cut, double high_pT_cut,
     double Psi_reaction_plane = 0.0);

  // Function that fills the v1, v2, and v3 profile histograms with data
  void fill_v1_v2_v3
    (double p_0, double p_x, double p_y, double p_z,
     // pass by reference to avoid copy; const to avoid changes
     const std::unique_ptr<TProfile>& v1_profile_histogram,
     const std::unique_ptr<TProfile>& v2_profile_histogram,
     const std::unique_ptr<TProfile>& v3_profile_histogram,
     double low_pT_cut, double high_pT_cut,
     double Psi_reaction_plane = 0.0);

  void fill_integrated_v1_v2_v3
    (double p_0, double p_x, double p_y, double p_z,
     // pass by reference to avoid copy; const to avoid changes
     const std::unique_ptr<TProfile>& v1_profile_histogram,
     const std::unique_ptr<TProfile>& v2_profile_histogram,
     const std::unique_ptr<TProfile>& v3_profile_histogram,
     double low_pT_cut, double high_pT_cut,
     double Psi_reaction_plane = 0.0);

  // Function that fits v1, v2, v3 profile histograms, creates plots, and saves fit data
  void fit_v1_v2_v3_make_plots_save_fit_data
    // pass by reference to avoid copy; const to avoid changes
    (const std::unique_ptr<TProfile>& p_v1,
     const std::unique_ptr<TProfile>& p_v2,
     const std::unique_ptr<TProfile>& p_v3,
     int n_events, char* impact_parameter_range_or_value, double time, Config cfg);

  // Function that saves contents of flow TProfiles into a data file (v1, v2, v3)
  void save_flow_vs_rapidity_data
    // pass by reference to avoid copy; const to avoid changes
    (const std::unique_ptr<TProfile>& p_v1_proton,
     const std::unique_ptr<TProfile>& p_v1_deuteron,
     const std::unique_ptr<TProfile>& p_v1_lambda,
     const std::unique_ptr<TProfile>& p_v1_pi_plus,
     const std::unique_ptr<TProfile>& p_v1_pi_minus,
     const std::unique_ptr<TProfile>& p_v1_kaon_plus,
     const std::unique_ptr<TProfile>& p_v1_kaon_minus,
     const std::unique_ptr<TProfile>& p_v2_proton,
     const std::unique_ptr<TProfile>& p_v2_deuteron,
     const std::unique_ptr<TProfile>& p_v2_lambda,
     const std::unique_ptr<TProfile>& p_v2_pi_plus,
     const std::unique_ptr<TProfile>& p_v2_pi_minus,
     const std::unique_ptr<TProfile>& p_v2_kaon_plus,
     const std::unique_ptr<TProfile>& p_v2_kaon_minus,
     const std::unique_ptr<TProfile>& p_v3_proton,
     const std::unique_ptr<TProfile>& p_v3_deuteron,
     const std::unique_ptr<TProfile>& p_v3_lambda,
     const std::unique_ptr<TProfile>& p_v3_pi_plus,
     const std::unique_ptr<TProfile>& p_v3_pi_minus,
     const std::unique_ptr<TProfile>& p_v3_kaon_plus,
     const std::unique_ptr<TProfile>& p_v3_kaon_minus,
     int n_events, char* impact_parameter_range_or_value, double time);

  // Function that saves total number of particle species (no cuts)
  void save_yield_data
  // pass by reference to avoid copy; const to avoid changes
  (const std::vector<int>& proton_count,
   const std::vector<int>& deuteron_count,
   const std::vector<int>& lambda_count,
   const std::vector<int>& pi_plus_count,
   const std::vector<int>& pi_minus_count,
   const std::vector<int>& kaon_plus_count,
   const std::vector<int>& kaon_minus_count,
   std::vector<double> current_t_steps,
   std::set<double> output_times,
   // passed only to identify the type of analysis:
   const std::unique_ptr<TProfile>& p_v1_proton,
   int n_events, int real_event_equivalent, Config cfg);

  // Function that creates plots of flow evolution and saves the data
  void plot_and_save_flow_time_evolution_data
    // pass by reference to avoid copy; const to avoid changes
    (const std::unique_ptr<TGraphErrors>& g_coll,
     const std::unique_ptr<TGraphErrors>& g_MF,
     const std::unique_ptr<TGraphErrors>& g_total,
     int n_events, int real_event_equivalent,
     double x_axis_min, double x_axis_max, double y_axis_min, double y_axis_max,
     int entry_step, Config cfg);

  // Function calculating flow from <cos(nphi)> using the ideal theoretical event plane
  void basic_flow
    (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo config_info,
     Config cfg, std::set<double> output_times = {});

  // Function extracting the mean-field and collision term contributions to flow binned in
  // rapidity
  void basic_flow_time_evolution_binned_in_y
    (const std::unique_ptr<ReadParticles>& ROOT_file, Config cfg);

  // Function extracting the mean-field and collision term contributions to flow
  // integrated over 4pi and at midrapidity
  void basic_flow_time_evolution_in_4pi_and_at_mid_y
    (const std::unique_ptr<ReadParticles>& ROOT_file, Config cfg);
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for private class members
  ////////////////////////////////////////////////////////////////////////////////////////

  // Only define return functions for members used outside of the class


  
 private:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Parameters of the simulation:
  // center-of-mass frame beam energy
  const double sqrts_;
  // fixed-target frame beam energy
  const double Ekin_;
  // center-of-mass frame beam rapidity
  const double y_beam_;
  ////////////////////////////////////////////////////////////////////////////////////////
  // Time evolution of flow:
  // whether to calculate flow only at the final time step
  const bool only_flow_at_final_output_;
  // number of time steps at which flow will be calculated
  const int number_of_time_steps_;  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Rapidity coverage:
  // whether to scale flow by the beam rapidity
  const bool scale_by_beam_rapidity_; 
  // number of bins along rapidity
  const int y_number_of_bins_;
  // considered rapidity range (used to build the TProfiles)
  const double y_min_;
  const double y_max_;
  // considered midrapidity range (used to calculate flow integrated at midrapidity)
  const double y_mid_min_;
  const double y_mid_max_;  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Transverse momentum pT cuts:
  const double proton_pT_min_;
  const double proton_pT_max_;
  const double deuteron_pT_min_;
  const double deuteron_pT_max_;
  const double lambda_pT_min_;
  const double lambda_pT_max_;
  const double pion_pT_min_;
  const double pion_pT_max_;
  const double kaon_pT_min_;
  const double kaon_pT_max_;
  ////////////////////////////////////////////////////////////////////////////////////////
  // Parameters for plotting:
  const bool user_axes_ranges_;
  double x_axis_lower_;
  double x_axis_upper_;
  double v1_y_axis_lower_;
  double v1_y_axis_upper_;
  double v2_y_axis_lower_;
  double v2_y_axis_upper_;
  double v3_y_axis_lower_;
  double v3_y_axis_upper_;
  double fitting_range_;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Observables:
  ///////////////////////////////////////////
  // Yields, binned in rapidity, for the entire range (4pi) and for midrapidity only:
  std::vector<std::unique_ptr<TH1D>> h_N_proton_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_deuteron_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_lambda_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_pi_plus_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_pi_minus_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_kaon_plus_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_kaon_minus_in_4pi_{};
  std::vector<std::unique_ptr<TH1D>> h_N_proton_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_deuteron_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_lambda_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_pi_plus_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_pi_minus_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_kaon_plus_at_mid_y_{};
  std::vector<std::unique_ptr<TH1D>> h_N_kaon_minus_at_mid_y_{};
  ///////////////////////////////////////////
  // Basic flow, binned in rapidity:
  // pT-integrated v1 from basic flow analysis
  std::vector<std::unique_ptr<TProfile>> p_v1_proton_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_deuteron_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_lambda_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_pi_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_pi_minus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_kaon_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_kaon_minus_basic_{};
  // pT-integrated v2 from basic flow analysis
  std::vector<std::unique_ptr<TProfile>> p_v2_proton_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_deuteron_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_lambda_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_pi_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_pi_minus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_kaon_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_kaon_minus_basic_{};
  // pT-integrated v3 from basic flow analysis
  std::vector<std::unique_ptr<TProfile>> p_v3_proton_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_deuteron_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_lambda_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_pi_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_pi_minus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_kaon_plus_basic_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_kaon_minus_basic_{};
  // pT- and y-integrated proton v1, v2, and v3 from basic flow analysis
  std::vector<std::unique_ptr<TProfile>> p_v1_proton_basic_integrated_in_4pi_{};
  std::vector<std::unique_ptr<TProfile>> p_v1_proton_basic_integrated_at_mid_y_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_proton_basic_integrated_in_4pi_{};
  std::vector<std::unique_ptr<TProfile>> p_v2_proton_basic_integrated_at_mid_y_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_proton_basic_integrated_in_4pi_{};
  std::vector<std::unique_ptr<TProfile>> p_v3_proton_basic_integrated_at_mid_y_{};
  // Track whether the basic flow analysis has been performed
  bool basic_flow_analysis_done_{}; // brace-initialized to false
  ///////////////////////////////////////////
  // Time evolution of basic flow, binned in rapidity:
  // time derivative of collision-term contribution only:d vn_{coll}/dt
  std::vector<std::unique_ptr<TGraphErrors>> g_dv1coll_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv2coll_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv3coll_dt_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_dv1coll_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv2coll_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv3coll_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv1coll_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv2coll_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv3coll_dt_proton_basic_integrated_at_mid_y_{};
  // time derivative of mean-field contribution only: dvn_{MF}/dt
  std::vector<std::unique_ptr<TGraphErrors>> g_dv1MF_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv2MF_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv3MF_dt_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_dv1MF_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv2MF_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv3MF_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv1MF_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv2MF_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv3MF_dt_proton_basic_integrated_at_mid_y_{};
  // time derivative of total flow: dvn_{total}/dt
  std::vector<std::unique_ptr<TGraphErrors>> g_dv1total_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv2total_dt_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_dv3total_dt_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_dv1total_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv2total_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv3total_dt_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_dv1total_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv2total_dt_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_dv3total_dt_proton_basic_integrated_at_mid_y_{};
  // integrated (from t=0) collision-term contribution only: vn_{coll}(t)
  std::vector<std::unique_ptr<TGraphErrors>> g_v1coll_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v2coll_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v3coll_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_v1coll_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v2coll_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v3coll_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v1coll_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v2coll_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v3coll_proton_basic_integrated_at_mid_y_{};
  // integrated (from t=0) mean-field contribution only: vn_{MF}(t)
  std::vector<std::unique_ptr<TGraphErrors>> g_v1MF_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v2MF_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v3MF_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_v1MF_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v2MF_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v3MF_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v1MF_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v2MF_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v3MF_proton_basic_integrated_at_mid_y_{};
  // integrated (from t=0) total flow: vn_{total}(t)
  std::vector<std::unique_ptr<TGraphErrors>> g_v1total_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v2total_proton_basic_{};
  std::vector<std::unique_ptr<TGraphErrors>> g_v3total_proton_basic_{};
  std::unique_ptr<TGraphErrors> g_v1total_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v2total_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v3total_proton_basic_integrated_in_4pi_{};
  std::unique_ptr<TGraphErrors> g_v1total_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v2total_proton_basic_integrated_at_mid_y_{};
  std::unique_ptr<TGraphErrors> g_v3total_proton_basic_integrated_at_mid_y_{};


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Fits:
  // v1_fits (basic flow)
  std::pair<double, double> v1_proton_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_deuteron_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_lambda_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_pi_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_pi_minus_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_kaon_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v1_kaon_minus_fit_ = {0.0, 0.0};
  // v2_fits (basic flow)
  std::pair<double, double> v2_proton_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_deuteron_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_lambda_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_pi_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_pi_minus_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_kaon_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v2_kaon_minus_fit_ = {0.0, 0.0};
  // v3_fits (basic flow)
  std::pair<double, double> v3_proton_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_deuteron_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_lambda_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_pi_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_pi_minus_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_kaon_plus_fit_ = {0.0, 0.0};
  std::pair<double, double> v3_kaon_minus_fit_ = {0.0, 0.0};
};

#endif
