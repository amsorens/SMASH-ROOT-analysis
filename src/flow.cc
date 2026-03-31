/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate collective flow.
//////////////////////////////////////////////////////////////////////////////////////////

#include "./flow.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>        // for std::unique_ptr
#include <sstream>
#include <stdio.h>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>

#include "./constants.h"
#include "./basic_functions.h"
#include "./plots.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"



void Flow::fill_N_vs_y
  (double p_0, double p_x, double p_y, double p_z,
   // pass by reference to avoid copy; const to avoid changes
   const std::unique_ptr<TH1D>& N_histogram,
   double low_pT_cut, double high_pT_cut) {
  
  const double pT = transverse_momentum (p_x, p_y);
  
  // Apply pT cut
  if ( (pT > low_pT_cut) && (pT < high_pT_cut) ) {
    // Get rapidity
    double y_r = rapidity(p_z, p_0);
    if ( scale_by_beam_rapidity_ ) {
      y_r = y_r/y_beam_;
    }

    // Fill the dN_dy histogram
    N_histogram->Fill(y_r);
  }
}



void Flow::fill_v1_v2_v3
  (double p_0, double p_x, double p_y, double p_z,
   // pass by reference to avoid copy; const to avoid changes
   const std::unique_ptr<TProfile>& v1_profile_histogram,
   const std::unique_ptr<TProfile>& v2_profile_histogram,
   const std::unique_ptr<TProfile>& v3_profile_histogram,
   double low_pT_cut, double high_pT_cut,
   double Psi_reaction_plane) {
  
  const double pT = transverse_momentum (p_x, p_y);
  // calculate in-situ to speed up:
  //const double pT = std::sqrt(p_x * p_x + p_y * p_y);
  
  // Apply pT cut
  if ( (pT > low_pT_cut) && (pT < high_pT_cut) ) {
    // Get rapidity
    double y_r = rapidity(p_z, p_0);
    if ( scale_by_beam_rapidity_ ) {
      y_r = y_r/y_beam_;
    }

    // Fill v1, v2, v3 profile histograms
    // Calculate phi w.r.t. the reaction plane
    const double phi_minus_reaction_plane = phi_angle(p_x, p_y, pT) - Psi_reaction_plane;
    v1_profile_histogram->Fill( y_r, cos(phi_minus_reaction_plane) );
    v2_profile_histogram->Fill( y_r, cos(2.0 * (phi_minus_reaction_plane)) );
    v3_profile_histogram->Fill( y_r, cos(3.0 * (phi_minus_reaction_plane)) );
  }
}



void Flow::fill_integrated_v1_v2_v3
  (double p_0, double p_x, double p_y, double p_z,
   // pass by reference to avoid copy; const to avoid changes
   const std::unique_ptr<TProfile>& v1_profile_histogram,
   const std::unique_ptr<TProfile>& v2_profile_histogram,
   const std::unique_ptr<TProfile>& v3_profile_histogram,
   double low_pT_cut, double high_pT_cut,
   double Psi_reaction_plane) {
  
  const double pT = transverse_momentum (p_x, p_y);
  // calculate in-situ to speed up:
  //const double pT = std::sqrt(p_x * p_x + p_y * p_y);
  
  // Apply pT cut
  if ( (pT > low_pT_cut) && (pT < high_pT_cut) ) {
    // Get rapidity
    double y_r = rapidity(p_z, p_0);
    if ( scale_by_beam_rapidity_ ) {
      y_r = y_r/y_beam_;
    }

    // Fill integrated v1, v2, v3 profile histograms
    // Calculate phi w.r.t. the reaction plane
    const double phi_minus_reaction_plane = phi_angle(p_x, p_y, pT) - Psi_reaction_plane;
    const double sgn_y = (y_r < 0.0) ? -1.0 : 1.0;
    v1_profile_histogram->Fill( y_r, sgn_y * cos(phi_minus_reaction_plane) );
    v2_profile_histogram->Fill( y_r, cos(2.0 * (phi_minus_reaction_plane)) );
    v3_profile_histogram->Fill( y_r, sgn_y * cos(3.0 * (phi_minus_reaction_plane)) );
  }
}



void Flow::fit_v1_v2_v3_make_plots_save_fit_data
  // pass by reference to avoid copy; const to avoid changes
  (const std::unique_ptr<TProfile>& p_v1,
   const std::unique_ptr<TProfile>& p_v2,
   const std::unique_ptr<TProfile>& p_v3,
   int n_events, char* impact_parameter_range_or_value, double time, Config cfg,
   const bool put_plots_in_separate_directory) {
  // If put_plots_in_separate_directory is true, then make sure the directory exists
  if (put_plots_in_separate_directory) {
    // Full path: Target_Directory + directory name
    std::filesystem::path folder_path =
      std::filesystem::path(Target_Directory) / flow_directory_name_;

    std::error_code ec;
    std::filesystem::create_directories(folder_path, ec);
    if (ec) {
      std::cerr << "Warning: could not create directory "
		<< folder_path.string() << ": " << ec.message() << "\n";
    }
  }
  
  // Establish whether the fits should be chatty
  const char* fit_options = cfg.verbose ? "R" : "RQ";
  
  // Establish the fitting range
  double fitting_range = 0.0;
  if ( scale_by_beam_rapidity_ ) {
    fitting_range = fitting_range_/y_beam_;
  } else {
    fitting_range = fitting_range_;
  }

  // Read off the particle name
  // (we can do this because all histogram names follow the same pattern)
  std::string name(p_v1->GetName());
  // Remove prefix "p_v1_" or "p_v2_"
  std::string prefix = "p_v1_";
  if (name.substr(0, prefix.size()) != prefix) {
    throw std::runtime_error("Histogram name does not start with expected prefix.");
  }
  std::string stripped = name.substr(prefix.size());
  // Remove trailing "_number" if present
  size_t last_underscore = stripped.find_last_of('_');
  if (last_underscore != std::string::npos &&
      last_underscore + 1 < stripped.size() &&
      std::all_of(stripped.begin() + last_underscore + 1, stripped.end(), ::isdigit)) {
    stripped = stripped.substr(0, last_underscore);
  }
  // Now split into particle name and analysis type by last underscore
  size_t split_position = stripped.find_last_of('_');
  if (split_position == std::string::npos) {
    throw std::runtime_error("Expected at least one underscore in stripped name.");
  }
  std::string particle_name = stripped.substr(0, split_position);
  std::string analysis_type = stripped.substr(split_position + 1);
  //std::cout << "particle_name = " << particle_name << std::endl;
  //std::cout << "analysis_type = " << analysis_type << std::endl;

  //
  // HACK!!!!
  //
  // TO DO: find an appropriate solution
  //
  if (particle_name == "proton_basic_integrated_in" || analysis_type == "4pi") {
    particle_name = "proton";
    analysis_type = "basic_integrated_in_4pi";
  } else if (particle_name == "proton_basic_integrated_at_mid_" || analysis_type == "y") {
    particle_name = "proton";
    analysis_type = "basic_integrated_at_mid_y";
  }
  //
  // HACK!!!!
  //
  //std::cout << "particle_name = " << particle_name << std::endl;
  //std::cout << "analysis_type = " << analysis_type << std::endl;

  

  // Establish particle_pmin, particle_pmax, particle_ymin, particle_ymax
  double particle_pT_min = 0.0;
  double particle_pT_max = 0.0;
  if ( particle_name == "proton" ) {
    particle_pT_min = proton_pT_min_;
    particle_pT_max = proton_pT_max_;
  } else if ( particle_name == "deuteron" ) {
    particle_pT_min = deuteron_pT_min_;
    particle_pT_max = deuteron_pT_max_;
  } else if ( particle_name == "lambda" ) {
    particle_pT_min = lambda_pT_min_;
    particle_pT_max = lambda_pT_max_;
  } else if ( particle_name == "pi_plus" ) {
    particle_pT_min = pion_pT_min_;
    particle_pT_max = pion_pT_max_;
  } else if ( particle_name == "pi_minus" ) {
    particle_pT_min = pion_pT_min_;
    particle_pT_max = pion_pT_max_;
  } else if ( particle_name == "kaon_plus" ) {
    particle_pT_min = kaon_pT_min_;
    particle_pT_max = kaon_pT_max_;
  } else if ( particle_name == "kaon_minus" ) {
    particle_pT_min = kaon_pT_min_;
    particle_pT_max = kaon_pT_max_;
  } 

  // Read off whether flow data is scaled by rapidity
  char scale_by_ybeam_info[Char_Array_Size];
  if ( scale_by_beam_rapidity_ ) {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity=true");
  } else {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity=false");
  }
  

  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Fit v1
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  // Set up the fitting function
  const char* fit_function_v1 = "[0] * x + [1] * x * x * x";
  TF1 *v1_fit = new TF1("v1_fit", fit_function_v1, -fitting_range, fitting_range);
  v1_fit->SetParameter(0, 0.05);
  v1_fit->SetParameter(1, 0.10);
  v1_fit->SetLineColor(6);
  v1_fit->SetLineStyle(2);
  if (cfg.verbose) {
    std::cout << "\nfit_v1_v2_v3_make_plots_save_fit_data:\n"
	      << "v1 fitting function: " << fit_function_v1 << std::endl;
  }
  // Perform the fit
  p_v1->Fit("v1_fit", fit_options);


  
  ///////////////////////////////////////////
  // Plot v1 and fit

  // Get the name for the plot and files
  char basic_file_name_v1[Char_Array_Size];
  snprintf(basic_file_name_v1, Char_Array_Size,
	   "%s%s%s_pTmin_%3.1f_pTmax_%3.1f_ymin_%3.2f_ymax_%3.2f_nEvents_%d_t_%.1f_fm_%s",
	   Target_Directory,
	   put_plots_in_separate_directory ? flow_directory_name_ : "",
	   p_v1->GetName(),
	   particle_pT_min, particle_pT_max, -fitting_range, fitting_range, n_events, time,
	   scale_by_ybeam_info);

  // Create canvas 
  TCanvas *canvas_v1 = new TCanvas("canvas_v1", basic_file_name_v1, 1200, 1200);
  set_canvas_properties(canvas_v1);

  // Graph v1
  p_v1->SetLineColor(2);
  p_v1->SetMarkerStyle(8);
  if (user_axes_ranges_) {
    const double x_axis_lower_aux = scale_by_beam_rapidity_ ?
      x_axis_lower_/y_beam_ : x_axis_lower_;
    const double x_axis_upper_aux = scale_by_beam_rapidity_ ?
      x_axis_upper_/y_beam_ : x_axis_upper_;
    p_v1->GetXaxis()->SetRangeUser(x_axis_lower_aux, x_axis_upper_aux);
    p_v1->GetYaxis()->SetRangeUser(v1_y_axis_lower_, v1_y_axis_upper_);
  }
  p_v1->Draw("Pz");
  // Graph the fit
  v1_fit->SetLineColor(4);
  v1_fit->Draw("same");
  
  // Create .C, .pdf, .png files
  create_canvas_files(canvas_v1, basic_file_name_v1);
  // Create a ROOT file
  create_a_histogram_ROOT_file(p_v1, basic_file_name_v1);
  // Delete canvas
  delete canvas_v1;



  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Fit v2
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  // Set up the fitting function
  const char* fit_function_v2 = "[0] + [1] * x * x";
  TF1 *v2_fit = new TF1("v2_fit", fit_function_v2, -fitting_range, fitting_range);
  v2_fit->SetParameter(0, -0.05);
  v2_fit->SetParameter(1, 0.10);
  v2_fit->SetLineColor(6);
  v2_fit->SetLineStyle(2);
  if (cfg.verbose) {
    std::cout << "\nfit_v1_v2_v3_make_plots_save_fit_data:\n"
	      << "v2 fitting function: " << fit_function_v2 << std::endl;
  }
  // Perform the fit
  p_v2->Fit("v2_fit", fit_options);  


  
  ///////////////////////////////////////////
  // Plot v2 and fit

  // Get the name for the plot and files
  char basic_file_name_v2[Char_Array_Size];
  snprintf(basic_file_name_v2, Char_Array_Size,
	   "%s%s%s_pTmin_%3.1f_pTmax_%3.1f_ymin_%3.2f_ymax_%3.2f_nEvents_%d_t_%.1f_fm_%s",
	   Target_Directory,
	   put_plots_in_separate_directory ? flow_directory_name_ : "",
	   p_v2->GetName(),
	   particle_pT_min, particle_pT_max, -fitting_range, fitting_range, n_events, time,
	   scale_by_ybeam_info);

  // Create canvas 
  TCanvas *canvas_v2 = new TCanvas("canvas_v2", basic_file_name_v2, 1200, 1200);
  set_canvas_properties(canvas_v2);

  // Graph v1
  p_v2->SetLineColor(2);
  p_v2->SetMarkerStyle(8);
  if (user_axes_ranges_) {
    const double x_axis_lower_aux = scale_by_beam_rapidity_ ?
      x_axis_lower_/y_beam_ : x_axis_lower_;
    const double x_axis_upper_aux = scale_by_beam_rapidity_ ?
      x_axis_upper_/y_beam_ : x_axis_upper_;
    p_v2->GetXaxis()->SetRangeUser(x_axis_lower_aux, x_axis_upper_aux);
    p_v2->GetYaxis()->SetRangeUser(v2_y_axis_lower_, v2_y_axis_upper_);
  }
  p_v2->Draw("Pz");
  // Graph the fit
  v2_fit->SetLineColor(4);
  v2_fit->Draw("same");

  // Create .C, .pdf, .png files
  create_canvas_files(canvas_v2, basic_file_name_v2);
  // Create a ROOT file
  create_a_histogram_ROOT_file(p_v2, basic_file_name_v2);
  // Delete canvas
  delete canvas_v2;



  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Fit v3
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  // Set up the fitting function
  const char* fit_function_v3 = "[0] * x + [1] * x * x * x";
  TF1 *v3_fit = new TF1("v3_fit", fit_function_v3, -fitting_range, fitting_range);
  v3_fit->SetParameter(0, 0.05);
  v3_fit->SetParameter(1, 0.10);
  v3_fit->SetLineColor(6);
  v3_fit->SetLineStyle(2);
  if (cfg.verbose) {
    std::cout << "\nfit_v1_v2_v3_make_plots_save_fit_data:\n"
	      << "v3 fitting function: " << fit_function_v3 << std::endl;
  }
  // Perform the fit
  p_v3->Fit("v3_fit", fit_options);


  
  ///////////////////////////////////////////
  // Plot v3 and fit

  // Get the name for the plot and files
  char basic_file_name_v3[Char_Array_Size];
  snprintf(basic_file_name_v3, Char_Array_Size,
	   "%s%s%s_pTmin_%3.1f_pTmax_%3.1f_ymin_%3.2f_ymax_%3.2f_nEvents_%d_t_%.1f_fm_%s",
	   Target_Directory,
	   put_plots_in_separate_directory ? flow_directory_name_ : "",
	   p_v3->GetName(),
	   particle_pT_min, particle_pT_max, -fitting_range, fitting_range, n_events, time,
	   scale_by_ybeam_info);

  // Create canvas 
  TCanvas *canvas_v3 = new TCanvas("canvas_v3", basic_file_name_v3, 1200, 1200);
  set_canvas_properties(canvas_v3);

  // Graph v3
  p_v3->SetLineColor(2);
  p_v3->SetMarkerStyle(8);
  if (user_axes_ranges_) {
    const double x_axis_lower_aux = scale_by_beam_rapidity_ ?
      x_axis_lower_/y_beam_ : x_axis_lower_;
    const double x_axis_upper_aux = scale_by_beam_rapidity_ ?
      x_axis_upper_/y_beam_ : x_axis_upper_;
    p_v3->GetXaxis()->SetRangeUser(x_axis_lower_aux, x_axis_upper_aux);
    p_v3->GetYaxis()->SetRangeUser(v3_y_axis_lower_, v3_y_axis_upper_);
  }
  p_v3->Draw("Pz");
  // Graph the fit
  v3_fit->SetLineColor(4);
  v3_fit->Draw("same");
  
  // Create .C, .pdf, .png files
  create_canvas_files(canvas_v3, basic_file_name_v3);
  // Create a ROOT file
  create_a_histogram_ROOT_file(p_v3, basic_file_name_v3);
  // Delete canvas
  delete canvas_v3;
  
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Update appropriate class members
  ////////////////////////////////////////////////////////////////////////////////////////

  if ( particle_name == "proton" ) {
    v1_proton_fit_.first = v1_fit->GetParameter(0);
    v1_proton_fit_.second = v1_fit->GetParError(0);
    v2_proton_fit_.first = v2_fit->GetParameter(0);
    v2_proton_fit_.second = v2_fit->GetParError(0);
    v3_proton_fit_.first = v3_fit->GetParameter(0);
    v3_proton_fit_.second = v3_fit->GetParError(0);
  } else if ( particle_name == "deuteron" ) {
    v1_deuteron_fit_.first = v1_fit->GetParameter(0);
    v1_deuteron_fit_.second = v1_fit->GetParError(0);
    v2_deuteron_fit_.first = v2_fit->GetParameter(0);
    v2_deuteron_fit_.second = v2_fit->GetParError(0);
    v3_deuteron_fit_.first = v3_fit->GetParameter(0);
    v3_deuteron_fit_.second = v3_fit->GetParError(0);
  } else if ( particle_name == "lambda" ) {
    v1_lambda_fit_.first = v1_fit->GetParameter(0);
    v1_lambda_fit_.second = v1_fit->GetParError(0);
    v2_lambda_fit_.first = v2_fit->GetParameter(0);
    v2_lambda_fit_.second = v2_fit->GetParError(0);
    v3_lambda_fit_.first = v3_fit->GetParameter(0);
    v3_lambda_fit_.second = v3_fit->GetParError(0);
  } else if ( particle_name == "pi_plus" ) {
    v1_pi_plus_fit_.first = v1_fit->GetParameter(0);
    v1_pi_plus_fit_.second = v1_fit->GetParError(0);
    v2_pi_plus_fit_.first = v2_fit->GetParameter(0);
    v2_pi_plus_fit_.second = v2_fit->GetParError(0);
    v3_pi_plus_fit_.first = v3_fit->GetParameter(0);
    v3_pi_plus_fit_.second = v3_fit->GetParError(0);
  } else if ( particle_name == "pi_minus" ) {
    v1_pi_minus_fit_.first = v1_fit->GetParameter(0);
    v1_pi_minus_fit_.second = v1_fit->GetParError(0);
    v2_pi_minus_fit_.first = v2_fit->GetParameter(0);
    v2_pi_minus_fit_.second = v2_fit->GetParError(0);
    v3_pi_minus_fit_.first = v3_fit->GetParameter(0);
    v3_pi_minus_fit_.second = v3_fit->GetParError(0);
  } 
  

  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Put the fit data into a local and into a global file
  ////////////////////////////////////////////////////////////////////////////////////////

  // Declare a string object with an appropriate name
  std::string Fit_data_file_name_basic(Target_Directory);
  Fit_data_file_name_basic += "00_fit_data_v1_v2_v3_";
  std::stringstream data_0;
  data_0 << particle_name;
  Fit_data_file_name_basic += data_0.str();
  Fit_data_file_name_basic += "_";
  std::stringstream data_1;
  data_1 << analysis_type;
  Fit_data_file_name_basic += data_1.str();


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Proceed for the local file option
  std::string Fit_data_local_file_name = Fit_data_file_name_basic;
  
  // Read off the pTcut
  if ( particle_name == "proton" ) {
    // these are protons
    if (proton_pT_min_ > 0.0) {
      Fit_data_local_file_name += "_pTmin=";
      std::stringstream data_pTcut_min;
      data_pTcut_min << proton_pT_min_;
      Fit_data_local_file_name += data_pTcut_min.str();
      Fit_data_local_file_name += "_pTmax=";
      std::stringstream data_pTcut_max;
      data_pTcut_max << proton_pT_max_;
      Fit_data_local_file_name += data_pTcut_max.str();
    } else {
      Fit_data_local_file_name += "_no_pTcut";
    }
  } else if ( particle_name == "deuteron" ) {
    // these are deuterons
    if (deuteron_pT_min_ > 0.0) {
      Fit_data_local_file_name += "_pTmin=";
      std::stringstream data_pTcut_min;
      data_pTcut_min << deuteron_pT_min_;
      Fit_data_local_file_name += data_pTcut_min.str();
      Fit_data_local_file_name += "_pTmax=";
      std::stringstream data_pTcut_max;
      data_pTcut_max << deuteron_pT_max_;
      Fit_data_local_file_name += data_pTcut_max.str();
    } else {
      Fit_data_local_file_name += "_no_pTcut";
    }
  } else if ( particle_name == "lambda" ) {
    // these are lambdas
    if (lambda_pT_min_ > 0.0) {
      Fit_data_local_file_name += "_pTmin=";
      std::stringstream data_pTcut_min;
      data_pTcut_min << lambda_pT_min_;
      Fit_data_local_file_name += data_pTcut_min.str();
      Fit_data_local_file_name += "_pTmax=";
      std::stringstream data_pTcut_max;
      data_pTcut_max << lambda_pT_max_;
      Fit_data_local_file_name += data_pTcut_max.str();
    } else {
      Fit_data_local_file_name += "_no_pTcut";
    }
  } else if ( (particle_name == "pi_plus") || (particle_name == "pi_minus") ) {
    // these are pions
    if (pion_pT_min_ > 0.0) {
      Fit_data_local_file_name += "_pTmin=";
      std::stringstream data_pTcut_min;
      data_pTcut_min << pion_pT_min_;
      Fit_data_local_file_name += data_pTcut_min.str();
      Fit_data_local_file_name += "_pTmax=";
      std::stringstream data_pTcut_max;
      data_pTcut_max << pion_pT_max_;
      Fit_data_local_file_name += data_pTcut_max.str();
    } else {
      Fit_data_local_file_name += "_no_pTcut";
    }
  } else if ( (particle_name == "kaon_plus") || (particle_name == "kaon_minus") ) {
    // these are kaons
    if (kaon_pT_min_ > 0.0) {
      Fit_data_local_file_name += "_pTmin=";
      std::stringstream data_pTcut_min;
      data_pTcut_min << kaon_pT_min_;
      Fit_data_local_file_name += data_pTcut_min.str();
      Fit_data_local_file_name += "_pTmax=";
      std::stringstream data_pTcut_max;
      data_pTcut_max << kaon_pT_max_;
      Fit_data_local_file_name += data_pTcut_max.str();
    } else {
      Fit_data_local_file_name += "_no_pTcut";
    }
  } else {
    std::cout << "\nfit_v1_v2_v3_make_plots_save_fit_data:\n"
	      << "this histogram is for particles for which no pT cut has been defined:"
	      << "\nparticle_name = " << particle_name
	      << "\nanalysis_type = " << analysis_type
	      << std::endl;
    std::cin.get();
  }

  // Read off fitted range
  Fit_data_local_file_name += "_ymin=";
  std::stringstream data_y_min;
  data_y_min << -fitting_range;
  Fit_data_local_file_name += data_y_min.str();
  Fit_data_local_file_name += "_ymax=";
  std::stringstream data_y_max;
  data_y_max << fitting_range;
  Fit_data_local_file_name += data_y_max.str();

  // Read off number of events
  Fit_data_local_file_name += "_nEvents=";
  std::stringstream data_nEvents;
  data_nEvents << n_events;
  Fit_data_local_file_name += data_nEvents.str();

  // Read off time
  Fit_data_local_file_name += "_t=";
  std::stringstream data_time;
  data_time << time;
  Fit_data_local_file_name += data_time.str();

  // Read off whether v1, v2, v3 are scaled by rapidity
  Fit_data_local_file_name += "_fm_scaled_rapidity=";
  if ( scale_by_beam_rapidity_ ) {
    Fit_data_local_file_name += "true";
  } else {
    Fit_data_local_file_name += "false";
  }

  // Add file extension
  Fit_data_local_file_name += ".txt";

  // Initialize the data file
  FILE *Fit_data_local;
  Fit_data_local = fopen(Fit_data_local_file_name.c_str(), "w");

  // Write relevant info into the file
  std::fprintf(Fit_data_local,
	       "#v1 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n",
	       v1_fit->GetChisquare(),
	       v1_fit->GetNDF(),
	       (v1_fit->GetChisquare())/(v1_fit->GetNDF()) );

  std::fprintf(Fit_data_local,
	       "#v2 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n",
	       v2_fit->GetChisquare(),
	       v2_fit->GetNDF(),
	       (v2_fit->GetChisquare())/(v2_fit->GetNDF()) );

  std::fprintf(Fit_data_local,
	       "#v3 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n#\n",
	       v3_fit->GetChisquare(),
	       v3_fit->GetNDF(),
	       (v3_fit->GetChisquare())/(v3_fit->GetNDF()) );

  std::fprintf(Fit_data_local,
	       "#sqrts        Ekin          y_beam         "
	       "dv1/dy(y=0)   +-            "
	       "v2(y=0)       +-            "
	       "dv3/dy(y=0)   +-            \n#\n");

  std::fprintf(Fit_data_local,
	       "%8.6f      %8.6f      %8.6f      %8.6f      %8.6f      "
	       "%8.6f      %8.6f      %8.6f      %8.6f      \n",
	       sqrts_, Ekin_, y_beam_,
	       v1_fit->GetParameter(0), v1_fit->GetParError(0),
	       v2_fit->GetParameter(0), v2_fit->GetParError(0),
	       v3_fit->GetParameter(0), v3_fit->GetParError(0) );
    

  fclose(Fit_data_local);



  ////////////////////////////////////////////////////////////////////////////////////////
  // Proceed for the global file option

  std::string Fit_data_global_file_name = "../../";
  Fit_data_global_file_name += Fit_data_file_name_basic;

  // Read off whether v1, v2, v3 are scaled by rapidity
  Fit_data_global_file_name += "_scaled_rapidity=";
  if ( scale_by_beam_rapidity_ ) {
    Fit_data_global_file_name += "true";
  } else {
    Fit_data_global_file_name += "false";
  }
  
  // Add file extension
  Fit_data_global_file_name += ".txt";

  // Initialize the data file
  FILE *Fit_data_global;

  // Check if the global file already exists
  std::ifstream file(Fit_data_global_file_name);
  bool file_exists = false;  
  if(!file)                  
    file_exists = false; // if the file was not found, then file=0, i.e. !file=1 or true.
  else
    file_exists = true; // If the file was found, then file is non-0.
  file.close();

  // Open the data file
  if ( !file_exists ) {
    // the file doesn't exist yet, create a new one and give it a header
    Fit_data_global = fopen(Fit_data_global_file_name.c_str(), "w");
    std::fprintf(Fit_data_global,
		 "# v1, v2, v3 fit data for various energies, centralities, and pT cuts\n"
		 "#\n"
		 "#sqrts         Ekin           y_beam         "
		 "dv1/dy(y=0)    +-             "
		 "v2(y=0)        +-             "
		 "dv3/dy(y=0)    +-            \n");
  } else {
    // the file exists, append the file
    Fit_data_global = fopen(Fit_data_global_file_name.c_str(), "a");
  }

  // Write relevant info into the file
  if ( particle_name == "proton" ) {
    std::fprintf(Fit_data_global,
		 "#\n#\n#pT_min = %4.2f, pT_max = %4.2f, impact parameter range = %s, "
		 "time = %.1f \n",
		 proton_pT_min_, proton_pT_max_, impact_parameter_range_or_value, time);
  } else if ( particle_name == "deuteron" ) {
    std::fprintf(Fit_data_global,
		 "#\n#\n#pT_min = %4.2f, pT_max = %4.2f, impact parameter range = %s, "
		 "time = %.1f \n",
		 deuteron_pT_min_, deuteron_pT_max_, impact_parameter_range_or_value, time);
  } else if ( particle_name == "lambda" ) {
    std::fprintf(Fit_data_global,
		 "#\n#\n#pT_min = %4.2f, pT_max = %4.2f, impact parameter range = %s, "
		 "time = %.1f \n",
		 lambda_pT_min_, lambda_pT_max_, impact_parameter_range_or_value, time);
  } else if ( (particle_name == "pi_plus") || (particle_name == "pi_minus") ) {
    std::fprintf(Fit_data_global,
		 "#\n#\n#pT_min = %4.2f, pT_max = %4.2f, impact parameter range = %s, "
		 "time = %.1f \n",
		 pion_pT_min_, pion_pT_max_, impact_parameter_range_or_value, time);
  }  else if ( (particle_name == "kaon_plus") || (particle_name == "kaon_minus") ) {
    std::fprintf(Fit_data_global,
		 "#\n#\n#pT_min = %4.2f, pT_max = %4.2f, impact parameter range = %s, "
		 "time = %.1f \n",
		 kaon_pT_min_, kaon_pT_max_, impact_parameter_range_or_value, time);
  }
    
  std::fprintf(Fit_data_global,
	       "#v1 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n",
	       v1_fit->GetChisquare(),
	       v1_fit->GetNDF(),
	       (v1_fit->GetChisquare())/(v1_fit->GetNDF()) );

  std::fprintf(Fit_data_global,
	       "#v2 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n",
	       v2_fit->GetChisquare(),
	       v2_fit->GetNDF(),
	       (v2_fit->GetChisquare())/(v2_fit->GetNDF()) );

  std::fprintf(Fit_data_global,
	       "#v3 fit: chi^2 = %8.6f    NDF = %d    chi^2/NDF = %8.6f \n",
	       v3_fit->GetChisquare(),
	       v3_fit->GetNDF(),
	       (v3_fit->GetChisquare())/(v3_fit->GetNDF()) );

  std::fprintf(Fit_data_global,
	       "%8.6f      %8.6f      %8.6f      %8.6f      %8.6f      "
	       "%8.6f      %8.6f      %8.6f      %8.6f      \n",
	       sqrts_, Ekin_, y_beam_,
	       v1_fit->GetParameter(0), v1_fit->GetParError(0),
	       v2_fit->GetParameter(0), v2_fit->GetParError(0),
	       v3_fit->GetParameter(0), v3_fit->GetParError(0) ); 

  fclose(Fit_data_global);


  
  // Delete fits
  delete v1_fit;
  delete v2_fit;
  delete v3_fit;
}



void Flow::save_flow_vs_rapidity_data
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
   int n_events, char* impact_parameter_range_or_value, double time) {
  ///////////////////////////////////////////
  // Establish the file name:
  // Identify the directory
  char Flow_data_directory[Char_Array_Size];
  snprintf(Flow_data_directory, Char_Array_Size, "%s", Target_Directory);
  
  // Identify the type of analysis (any p_v1 will do, we use proton)
  std::string name(p_v1_proton->GetName());
  // Step 1: Remove trailing "_number" if it exists
  size_t last_underscore = name.find_last_of('_');
  bool has_number_suffix =
    (last_underscore != std::string::npos &&
     last_underscore + 1 < name.size() &&
     std::all_of(name.begin() + last_underscore + 1, name.end(),::isdigit));
  std::string trimmed = has_number_suffix ? name.substr(0, last_underscore) : name;
  // Step 2: Extract word after the last underscore (this is the analysis type)
  size_t analysis_underscore = trimmed.find_last_of('_');
  if (analysis_underscore == std::string::npos) {
    throw std::runtime_error("Expected at least one underscore in histogram name.");
  }
  std::string analysis_type = trimmed.substr(analysis_underscore + 1);
  
  // Read off whether flow data is scaled by rapidity
  char scale_by_ybeam_info[Char_Array_Size];
  if ( scale_by_beam_rapidity_ ) {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity=true");
  } else {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity=false");
  }
  // Put all data file name info together
  char Flow_data_full_file_name[Char_Array_Size];
  snprintf(Flow_data_full_file_name, Char_Array_Size,
	   "%s00_flow_data_v1_v2_v3_%s_sqrts=%.2f_nEvents=%d_t=%.1f_fm_%s.txt",
	   Flow_data_directory, analysis_type.c_str(), sqrts_, n_events, time,
	   scale_by_ybeam_info);
  FILE *Flow_data;
  Flow_data = fopen(Flow_data_full_file_name, "w");

  // Add header
  std::fprintf(Flow_data, "# Flow data                \n#\n"
	       "# impact parameter range = %s fm \n"
	       "#   proton pTmin = %3.2f, pTmax = %3.2f \n"
	       "# deuteron pTmin = %3.2f, pTmax = %3.2f \n"
	       "#   Lambda pTmin = %3.2f, pTmax = %3.2f \n"
	       "#     pion pTmin = %3.2f, pTmax = %3.2f \n"
	       "#     kaon pTmin = %3.2f, pTmax = %3.2f \n#\n",
	       impact_parameter_range_or_value,
	       proton_pT_min_, proton_pT_max_, deuteron_pT_min_, deuteron_pT_max_,
	       lambda_pT_min_, lambda_pT_max_, pion_pT_min_, pion_pT_max_,
	       kaon_pT_min_, kaon_pT_max_);

  // Add column labels
  std::fprintf(Flow_data, "#y           "
	       "v1 p        +-          v2 p        +-          v3 p        +-          "
	       "v1 deuteron +-          v2 deuteron +-          v3 deuteron +-          "
	       "v1 Lambda   +-          v2 Lambda   +-          v3 Lambda   +-          "
	       "v1 pi plus  +-          v2 pi plus  +-          v3 pi plus  +-          "
	       "v1 pi minus +-          v2 pi minus +-          v3 pi minus +-          "
	       "v1 K+       +-          v2 K+       +-          v3 K+       +-          "
	       "v1 K-       +-          v2 K-       +-          v3 K-       +-          "
	       "\n#\n");


  
  ///////////////////////////////////////////
  // Add data
  // note: TProfile indices start at 1
  for (int i = 1; i < (y_number_of_bins_ + 1); i++){
    // the bin number is the same for all cases
    const int bin_number_proton_v1 = p_v1_proton->GetBin(i);
    // the rapidity is the same for all cases
    const double y_r = p_v1_proton->GetBinCenter(bin_number_proton_v1);

    std::pair<double, double> v1_proton = bin_value_and_error(p_v1_proton, i);
    std::pair<double, double> v2_proton = bin_value_and_error(p_v2_proton, i);
    std::pair<double, double> v3_proton = bin_value_and_error(p_v3_proton, i);
    std::pair<double, double> v1_deuteron = bin_value_and_error(p_v1_deuteron, i);
    std::pair<double, double> v2_deuteron = bin_value_and_error(p_v2_deuteron, i);
    std::pair<double, double> v3_deuteron = bin_value_and_error(p_v3_deuteron, i);
    std::pair<double, double> v1_lambda = bin_value_and_error(p_v1_lambda, i);
    std::pair<double, double> v2_lambda = bin_value_and_error(p_v2_lambda, i);
    std::pair<double, double> v3_lambda = bin_value_and_error(p_v3_lambda, i);
    std::pair<double, double> v1_pi_plus = bin_value_and_error(p_v1_pi_plus, i);
    std::pair<double, double> v2_pi_plus = bin_value_and_error(p_v2_pi_plus, i);
    std::pair<double, double> v3_pi_plus = bin_value_and_error(p_v3_pi_plus, i);
    std::pair<double, double> v1_pi_minus = bin_value_and_error(p_v1_pi_minus, i);
    std::pair<double, double> v2_pi_minus = bin_value_and_error(p_v2_pi_minus, i);
    std::pair<double, double> v3_pi_minus = bin_value_and_error(p_v3_pi_minus, i);
    std::pair<double, double> v1_K_plus = bin_value_and_error(p_v1_kaon_plus, i);
    std::pair<double, double> v2_K_plus = bin_value_and_error(p_v2_kaon_plus, i);
    std::pair<double, double> v3_K_plus = bin_value_and_error(p_v3_kaon_plus, i);
    std::pair<double, double> v1_K_minus = bin_value_and_error(p_v1_kaon_minus, i);
    std::pair<double, double> v2_K_minus = bin_value_and_error(p_v2_kaon_minus, i);
    std::pair<double, double> v3_K_minus = bin_value_and_error(p_v3_kaon_minus, i);

    std::fprintf(Flow_data, "%8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    \n",
		 y_r,
		 v1_proton.first, v1_proton.second,
		 v2_proton.first, v2_proton.second,
		 v3_proton.first, v3_proton.second,
		 v1_deuteron.first, v1_deuteron.second,
		 v2_deuteron.first, v2_deuteron.second,
		 v3_deuteron.first, v3_deuteron.second,
		 v1_lambda.first, v1_lambda.second,
		 v2_lambda.first, v2_lambda.second,
		 v3_lambda.first, v3_lambda.second,
		 v1_pi_plus.first, v1_pi_plus.second,
		 v2_pi_plus.first, v2_pi_plus.second,
		 v3_pi_plus.first, v3_pi_plus.second,
		 v1_pi_minus.first, v1_pi_minus.second,
		 v2_pi_minus.first, v2_pi_minus.second,
		 v3_pi_minus.first, v3_pi_minus.second,
		 v1_K_plus.first, v1_K_plus.second,
		 v2_K_plus.first, v2_K_plus.second,
		 v3_K_plus.first, v3_K_plus.second,
		 v1_K_minus.first, v1_K_minus.second,
		 v2_K_minus.first, v2_K_minus.second,
		 v3_K_minus.first, v3_K_minus.second);    
  }

  fclose(Flow_data);
}



void Flow::save_yield_data
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
   int n_events, int real_event_equivalent, Config cfg) {
  ///////////////////////////////////////////
  // Establish the file name:
  // Identify the directory
  char Yield_data_directory[Char_Array_Size];
  snprintf(Yield_data_directory, Char_Array_Size, "%s", Target_Directory);

  // Identify the type of analysis (any p_v1 will do, we use proton)
  std::string name(p_v1_proton->GetName());
  // Step 1: Remove trailing "_number" if it exists
  size_t last_underscore = name.find_last_of('_');
  bool has_number_suffix =
    (last_underscore != std::string::npos &&
     last_underscore + 1 < name.size() &&
     std::all_of(name.begin() + last_underscore + 1, name.end(),::isdigit));
  std::string trimmed = has_number_suffix ? name.substr(0, last_underscore) : name;
  // Step 2: Extract word after the last underscore (this is the analysis type)
  size_t analysis_underscore = trimmed.find_last_of('_');
  if (analysis_underscore == std::string::npos) {
    throw std::runtime_error("Expected at least one underscore in histogram name.");
  }
  std::string analysis_type = trimmed.substr(analysis_underscore + 1);
  
  // Put all data file name info together
  char Yield_data_full_file_name[Char_Array_Size];
  snprintf(Yield_data_full_file_name, Char_Array_Size,
	   "%s00_yield_data_from_%s_flow_analysis_sqrts=%.2f_nEvents=%d.txt",
	   Yield_data_directory, analysis_type.c_str(), sqrts_, n_events);
  FILE *Yield_data;
  Yield_data = fopen(Yield_data_full_file_name, "w");
  
  // Add header
  std::fprintf(Yield_data, "# Yield data (per equivalent real event) from %d events,\n"
	       "# equivalent to %d real events\n"
	       "# (no cuts) \n#\n",
	       n_events, real_event_equivalent);

  // Add column labels
  std::fprintf(Yield_data, "#t[fm]       "
	       "proton      deuteron    Lambda      "
	       "pi plus     pi minus    K+          K-          \n#\n");


  
  ///////////////////////////////////////////
  // Print and save data
  for (int i = 0; i < number_of_time_steps_; i++) {
    if (output_times.size() > 0) {
      // Print and save data only for the selected times:
      if ( output_times.count( current_t_steps[i] ) ) {
	// Evaluate particle counts per real event equivalent (REE)
	const double proton_count_REE = 1.0*proton_count[i]/real_event_equivalent;
	const double deuteron_count_REE = 1.0*deuteron_count[i]/real_event_equivalent;
	const double lambda_count_REE = 1.0*lambda_count[i]/real_event_equivalent;
	const double pi_plus_count_REE = 1.0*pi_plus_count[i]/real_event_equivalent;
	const double pi_minus_count_REE = 1.0*pi_minus_count[i]/real_event_equivalent;
	const double kaon_plus_count_REE = 1.0*kaon_plus_count[i]/real_event_equivalent;
	const double kaon_minus_count_REE = 1.0*kaon_minus_count[i]/real_event_equivalent;
	if (cfg.verbose) {
	  std::cout << "\n"
		    << "\n            t = " << current_t_steps[i]
		    << "\n"
		    << "\n     proton # = " << proton_count[i]
		    << "\n   deuteron # = " << deuteron_count[i]
		    << "\n     lambda # = " << lambda_count[i]
		    << "\n    pi plus # = " << pi_plus_count[i]
		    << "\n   pi minus # = " << pi_minus_count[i]
		    << "\n  kaon plus # = " << kaon_plus_count[i]
		    << "\n kaon minus # = " << kaon_minus_count[i]
		    << "\n"
		    << "\n    <proton #>/event = " << proton_count_REE
		    << "\n  <deuteron #>/event = " << deuteron_count_REE
		    << "\n    <lambda #>/event = " << lambda_count_REE
		    << "\n   <pi plus #>/event = " << pi_plus_count_REE
		    << "\n  <pi minus #>/event = " << pi_minus_count_REE
		    << "\n <kaon plus #>/event = " << kaon_plus_count_REE
		    << "\n<kaon minus #>/event = " << kaon_minus_count_REE
		    << "\n" << std::endl;
	}

	std::fprintf(Yield_data, "%6.1f     "
		     "%8.2f  %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    \n",
		     current_t_steps[i],
		     proton_count_REE, deuteron_count_REE, lambda_count_REE,
		     pi_plus_count_REE, pi_minus_count_REE,
		     kaon_plus_count_REE, kaon_minus_count_REE); 
      }
    } else {
      // Print and save data for all possible output times:
      // Evaluate particle counts per real event equivalent
      const double proton_count_REE = 1.0*proton_count[i]/real_event_equivalent;
      const double deuteron_count_REE = 1.0*deuteron_count[i]/real_event_equivalent;
      const double lambda_count_REE = 1.0*lambda_count[i]/real_event_equivalent;
      const double pi_plus_count_REE = 1.0*pi_plus_count[i]/real_event_equivalent;
      const double pi_minus_count_REE = 1.0*pi_minus_count[i]/real_event_equivalent;
      const double kaon_plus_count_REE = 1.0*kaon_plus_count[i]/real_event_equivalent;
      const double kaon_minus_count_REE = 1.0*kaon_minus_count[i]/real_event_equivalent;
      if (cfg.verbose) {
	std::cout << "\n"
		  << "\n            t = " << current_t_steps[i]
		  << "\n"
		  << "\n     proton # = " << proton_count[i]
		  << "\n   deuteron # = " << deuteron_count[i]
		  << "\n     lambda # = " << lambda_count[i]
		  << "\n    pi plus # = " << pi_plus_count[i]
		  << "\n   pi minus # = " << pi_minus_count[i]
		  << "\n  kaon plus # = " << kaon_plus_count[i]
		  << "\n kaon minus # = " << kaon_minus_count[i]
		  << "\n"
		  << "\n    <proton #>/event = " << proton_count_REE
		  << "\n  <deuteron #>/event = " << deuteron_count_REE
		  << "\n    <lambda #>/event = " << lambda_count_REE
		  << "\n   <pi plus #>/event = " << pi_plus_count_REE
		  << "\n  <pi minus #>/event = " << pi_minus_count_REE
		  << "\n <kaon plus #>/event = " << kaon_plus_count_REE
		  << "\n<kaon minus #>/event = " << kaon_minus_count_REE
		  << "\n" << std::endl;
      }

      std::fprintf(Yield_data, "%6.1f      "
		   "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    \n",
		   current_t_steps[i],
		   proton_count_REE, deuteron_count_REE, lambda_count_REE,
		   pi_plus_count_REE, pi_minus_count_REE,
		   kaon_plus_count_REE, kaon_minus_count_REE);  
    }
  }

  fclose(Yield_data);
}



void Flow::plot_and_save_flow_time_evolution_data
  // pass by reference to avoid copy; const to avoid changes
  (const std::unique_ptr<TGraphErrors>& g_coll,
   const std::unique_ptr<TGraphErrors>& g_MF,
   const std::unique_ptr<TGraphErrors>& g_total,
   int n_events, int real_event_equivalent,
   double x_axis_min, double x_axis_max, double y_axis_min, double y_axis_max,
   int entry_step, Config cfg, const bool put_plots_in_separate_directory) {
  // If put_plots_in_separate_directory is true, then make sure the directory exists
  if (put_plots_in_separate_directory) {
    // Full path: Target_Directory + directory name
    std::filesystem::path folder_path =
      std::filesystem::path(Target_Directory) / flow_evolution_directory_name_;

    std::error_code ec;
    std::filesystem::create_directories(folder_path, ec);
    if (ec) {
      std::cerr << "Warning: could not create directory "
		<< folder_path.string() << ": " << ec.message() << "\n";
    }
  }
  
  ///////////////////////////////////////////
  // Establish the file name:
  // Identify the directory
  char Flow_evolution_data_directory[Char_Array_Size];
  snprintf(Flow_evolution_data_directory, Char_Array_Size, "%s%s", Target_Directory,
	   put_plots_in_separate_directory ? flow_evolution_directory_name_ : "");
  
  // Identify the type of analysis, name structure: "g_dv1colldt_proton_basic_y_%.3f"
  // (any g_vn will do, we use g_coll)
  std::string full_name(g_coll->GetName());
  // Step 1: Remove the "g_" prefix
  if (full_name.substr(0, 2) != "g_") {
    throw std::runtime_error("Graph name must start with 'g_'.");
  }
  std::string name = full_name.substr(2);
  // Step 2: Optionally extract "_y_<number>"
  /*
  std::string ybin;
  size_t y_pos = name.rfind("_y_");
  if (y_pos != std::string::npos) {
    ybin = name.substr(y_pos + 3);  // part after "_y_"
    name = name.substr(0, y_pos);   // part before "_y_"
  }
  */
  std::string ybin;
  size_t y_pos = name.rfind("_y_");
  if (y_pos != std::string::npos) {
    std::string after_y = name.substr(y_pos + 3);
    // Check that what's after "_y_" looks like a number (optional sign, digits);
    // otherwise at may be the end of "at_mid_y"
    if ((!after_y.empty() && std::isdigit(after_y[0]))
	|| after_y[0] == '-' || after_y[0] == '+') {
      ybin = after_y;
      name = name.substr(0, y_pos);
    }
  }
  // Step 3: Split remaining name by underscores
  std::vector<std::string> tokens;
  std::stringstream ss(name);
  std::string item;
  while (std::getline(ss, item, '_')) {
    if (!item.empty()) tokens.push_back(item);
  }
  if (tokens.size() < 3) {
    throw std::runtime_error("Name must have at least one observable, particle, and "
			     "analysis fields.");
  }
  // Step 4: Extract last two fields
  std::string analysis_type = tokens.back();
  tokens.pop_back();
  std::string particle_name = tokens.back();
  tokens.pop_back();
  // Step 5: Combine the remaining fields as observable
  std::string observable;
  for (size_t i = 0; i < tokens.size(); ++i) {
    observable += tokens[i];
    if (i + 1 < tokens.size()) observable += "_";
  }
  // Step 6: Clean suffix tags in observable
  for (const std::string& tag : {"coll", "MF", "total"}) {
    size_t pos = observable.find(tag);
    if (pos != std::string::npos) {
      observable.erase(pos, tag.length());
      break;
    }
  }

  

  // Read off the pTcut
  char pT_cut_info[Char_Array_Size];
  if ( particle_name == "proton" ) {
    // these are protons
    if (proton_pT_min_ > 0.0) {
      const std::string pT_min = double_to_string_with_the_word_point(proton_pT_min_, 1);
      const std::string pT_max = double_to_string_with_the_word_point(proton_pT_max_, 1);
      snprintf(pT_cut_info, Char_Array_Size, "pTmin_%s_pTmax_%s",
	       pT_min.c_str(), pT_max.c_str());
    } else {
      snprintf(pT_cut_info, Char_Array_Size, "no_pTcut");
    }
  } else if ( particle_name == "deuteron" ) {
    // these are deuterons
    if (deuteron_pT_min_ > 0.0) {
      const std::string pT_min =
	double_to_string_with_the_word_point(deuteron_pT_min_, 1);
      const std::string pT_max =
	double_to_string_with_the_word_point(deuteron_pT_max_, 1);
      snprintf(pT_cut_info, Char_Array_Size, "pTmin_%s_pTmax_%s",
	       pT_min.c_str(), pT_max.c_str());
    } else {
      snprintf(pT_cut_info, Char_Array_Size, "no_pTcut");
    }
  } else if ( particle_name == "lambda" ) {
    // these are lambdas
    if (lambda_pT_min_ > 0.0) {
      const std::string pT_min = double_to_string_with_the_word_point(lambda_pT_min_, 1);
      const std::string pT_max = double_to_string_with_the_word_point(lambda_pT_max_, 1);
      snprintf(pT_cut_info, Char_Array_Size, "pTmin_%s_pTmax_%s",
	       pT_min.c_str(), pT_max.c_str());
    } else {
      snprintf(pT_cut_info, Char_Array_Size, "no_pTcut");
    }
  } else if ( (particle_name == "pi_plus") || (particle_name == "pi_minus") ) {
    // these are pions
    if (pion_pT_min_ > 0.0) {
      const std::string pT_min = double_to_string_with_the_word_point(pion_pT_min_, 1);
      const std::string pT_max = double_to_string_with_the_word_point(pion_pT_max_, 1);
      snprintf(pT_cut_info, Char_Array_Size, "pTmin_%s_pTmax_%s",
	       pT_min.c_str(), pT_max.c_str());
    } else {
      snprintf(pT_cut_info, Char_Array_Size, "no_pTcut");
    }
  } else if ( (particle_name == "kaon_plus") || (particle_name == "kaon_minus") ) {
    // these are kaons
    if (kaon_pT_min_ > 0.0) {
      const std::string pT_min = double_to_string_with_the_word_point(kaon_pT_min_, 1);
      const std::string pT_max = double_to_string_with_the_word_point(kaon_pT_max_, 1);
      snprintf(pT_cut_info, Char_Array_Size, "pTmin_%s_pTmax_%s",
	       pT_min.c_str(), pT_max.c_str());
    } else {
      snprintf(pT_cut_info, Char_Array_Size, "no_pTcut");
    }
  } else {
    std::cout << "\nplot_and_save_flow_time_evolution_data:\n"
	      << "this tgraph is for particles for which no pT cut has been defined"
	      << std::endl;
    std::cin.get();
  }
 
  // Read off whether flow data is scaled by rapidity
  char scale_by_ybeam_info[Char_Array_Size];
  if ( scale_by_beam_rapidity_ ) {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity_true");
  } else {
    snprintf(scale_by_ybeam_info, Char_Array_Size, "scaled_rapidity_false");
  }

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Create auxiliary TGraphs to plot a limited number of points:

  ///////////////////////////////////////////
  // coll
  char g_coll_aux_name[Char_Array_Size];
  snprintf(g_coll_aux_name, Char_Array_Size, "%s_each_%dth_step",
	   g_coll->GetName(), entry_step);
  auto g_coll_aux = std::make_unique<TGraphErrors>();
  g_coll_aux->SetName(g_coll_aux_name);
  g_coll_aux->SetTitle(g_coll_aux_name);
  // Copy points
  const int n_coll = g_coll->GetN();
  for (int i = 0; i < n_coll; i += entry_step) {
    double x = 0.0;
    double y = 0.0;
    g_coll->GetPoint(i, x, y);
    const int new_index = i / entry_step;
    g_coll_aux->SetPoint(new_index, x, y);
    g_coll_aux->SetPointError(new_index, g_coll->GetErrorX(i), g_coll->GetErrorY(i));
  }
  // MF
  char g_MF_aux_name[Char_Array_Size];
  snprintf(g_MF_aux_name, Char_Array_Size, "%s_each_%dth_step",
	   g_MF->GetName(), entry_step);
  auto g_MF_aux = std::make_unique<TGraphErrors>();
  g_MF_aux->SetName(g_MF_aux_name);
  g_MF_aux->SetTitle(g_MF_aux_name);
  // Copy points
  const int n_MF = g_MF->GetN();
  const int offset_MF = 1;
  for (int i = offset_MF; i < n_MF; i += entry_step) {
    double x = 0.0;
    double y = 0.0;
    g_MF->GetPoint(i, x, y);
    const int new_index = (i - offset_MF) / entry_step;
    g_MF_aux->SetPoint(new_index, x, y);
    g_MF_aux->SetPointError(new_index, g_MF->GetErrorX(i), g_MF->GetErrorY(i));
  }
  // total
  char g_total_aux_name[Char_Array_Size];
  snprintf(g_total_aux_name, Char_Array_Size, "%s_each_%dth_step",
	   g_total->GetName(), entry_step);
  auto g_total_aux = std::make_unique<TGraphErrors>();
  g_total_aux->SetName(g_total_aux_name);
  g_total_aux->SetTitle(g_total_aux_name);
  // Copy points
  const int n_total = g_total->GetN();
  const int offset_total = 2;
  for (int i = offset_total; i < n_total; i += entry_step) {
    double x = 0.0;
    double y = 0.0;
    g_total->GetPoint(i, x, y);
    const int new_index = (i - offset_total) / entry_step;
    g_total_aux->SetPoint(new_index, x, y);
    g_total_aux->SetPointError(new_index, g_total->GetErrorX(i), g_total->GetErrorY(i));
  }

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Plot TGraphs

  // Define auxiliary variables
  const std::string sqrts_string = double_to_string_with_the_word_point(sqrts_, 1);
  const std::string ybin_string =
    !ybin.empty() ? double_to_string_with_the_word_point(std::stod(ybin), 3) : "";
  const std::string ybin_info_tgraph = !ybin.empty() ? ybin_string : "";
  
  // Define tgraph file name
  char Flow_evolution_tgraph_file_name[Char_Array_Size];
  snprintf(Flow_evolution_tgraph_file_name, Char_Array_Size,
	   "%sg_%s_%s_%s_sqrts_%s_%s_nEvents_%d_%s%s",
	   Flow_evolution_data_directory, observable.c_str(), particle_name.c_str(),
	   analysis_type.c_str(), sqrts_string.c_str(), pT_cut_info, n_events,
	   ybin_info_tgraph.c_str(), scale_by_ybeam_info);

  // Create legend entries and an auxiliary vector
  std::vector<const char*> legend_entries = {"coll", "MF", "total"};
  std::vector<TGraphErrors*>
    tgraph_vector = {g_coll_aux.get(), g_MF_aux.get(), g_total_aux.get()};

  // Plot and save tgraphs
  plot_and_save_TGraph(tgraph_vector,
		       Flow_evolution_tgraph_file_name,
		       "t [fm/c]",
		       observable.c_str(),
		       legend_entries,
		       false, false,
		       true, x_axis_min, x_axis_max, y_axis_min, y_axis_max);
  


  ////////////////////////////////////////////////////////////////////////////////////////
  // Save data from TGraphs

  // Define auxiliary variable
  std::string ybin_info_txt = !ybin.empty() ? ("y=" + ybin_string + "_") : "";
  // Put all data file name info together
  char Flow_evolution_data_full_file_name[Char_Array_Size];
  snprintf(Flow_evolution_data_full_file_name, Char_Array_Size,
	   "%s00_flow_evolution_data_%s_%s_%s_sqrts=%s_%s_nEvents=%d_%s%s.txt",
	   Flow_evolution_data_directory, observable.c_str(), particle_name.c_str(),
	   analysis_type.c_str(), sqrts_string.c_str(), pT_cut_info, n_events,
	   ybin_info_txt.c_str(), scale_by_ybeam_info);
  //std::cin.get();
  FILE *Flow_evolution_data;
  Flow_evolution_data = fopen(Flow_evolution_data_full_file_name, "w");

  // Add header
  std::fprintf(Flow_evolution_data, "# Flow evolution data for %s \n"
	       "# from %d events, equivalent to %d real events\n#\n",
	       particle_name.c_str(), n_events, real_event_equivalent);

  // Add column labels
  std::fprintf(Flow_evolution_data, "#t           "
	       "%5s coll   +-      %5s MF     +-       %5s total      +-      \n#\n",
	       observable.c_str(), observable.c_str(), observable.c_str());
  
  ///////////////////////////////////////////
  // Add data
  for (int i = 0; i < g_coll->GetN(); i++){
    // the position on the x axis is the same for all tgraphs
    const double current_t = g_coll->GetX()[i];
    const double y_coll = g_coll->GetY()[i];
    const double y_coll_err = g_coll->GetEY()[i];
    const double y_MF = g_MF->GetY()[i];
    const double y_MF_err = g_MF->GetEY()[i];
    const double y_total = g_total->GetY()[i];
    const double y_total_err = g_total->GetEY()[i];

    std::fprintf(Flow_evolution_data, "%8.6f    "
		 "%8.6f    %8.6f    %8.6f    %8.6f    %8.6f    %8.6f    \n",
		 current_t,
		 y_coll, y_coll_err, y_MF, y_MF_err, y_total, y_total_err);   
  }

  fclose(Flow_evolution_data);
}



void Flow::basic_flow
  (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo config_info,
   Config cfg, std::set<double> output_times)
{
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform basic flow analysis, i.e., fill profile histograms with cos(nphi)
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Starting basic flow analysis                        :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n"
	    << color::RESET << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Enable only needed branches for this analysis:
  // Disable all
  ROOT_file->fChain->SetBranchStatus("*", 0);
  // Enable needed ones:
  ROOT_file->fChain->SetBranchStatus("t", 1);
  ROOT_file->fChain->SetBranchStatus("current_t", 1);
  ROOT_file->fChain->SetBranchStatus("npart", 1);
  ROOT_file->fChain->SetBranchStatus("pdgcode", 1);
  ROOT_file->fChain->SetBranchStatus("p0", 1);
  ROOT_file->fChain->SetBranchStatus("px", 1);
  ROOT_file->fChain->SetBranchStatus("py", 1);
  ROOT_file->fChain->SetBranchStatus("pz", 1);
  

  
  ///////////////////////////////////////////
  // Count particles
  std::vector<int> proton_count(number_of_time_steps_, 0);
  std::vector<int> deuteron_count(number_of_time_steps_, 0);
  std::vector<int> lambda_count(number_of_time_steps_, 0);
  std::vector<int> pi_plus_count(number_of_time_steps_, 0);
  std::vector<int> pi_minus_count(number_of_time_steps_, 0);
  std::vector<int> kaon_plus_count(number_of_time_steps_, 0);
  std::vector<int> kaon_minus_count(number_of_time_steps_, 0);


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Loop over all used entries

  const int progress_message_threshold = 1000;
  std::cout << "Print progress message every "
	    << progress_message_threshold << " entries:" << std::endl;

  const int entries_per_event =
    ROOT_file->n_event_time_steps() * ROOT_file->n_ensembles();
  
  for (int i_tree_entry = 0; i_tree_entry < ROOT_file->n_entries(); i_tree_entry++) {
    // Print out a message every X entries
    if ( (i_tree_entry % progress_message_threshold) == 0 ) {
      std::cout << "basic flow analysis: loading i_tree_entry = "
		<< i_tree_entry << std::endl;
    }

    // Load the TTree data
    Long64_t Long64_t_for_entry = ROOT_file->LoadTree(i_tree_entry);
    if (Long64_t_for_entry < 0) {
      std::cout << "i_tree_entry = " << i_tree_entry << std::endl;
      throw std::runtime_error("Failed to load the TTree at the above entry.");
    }
    long long int l_l_int_for_entry = ROOT_file->GetEntry(i_tree_entry);

    

    //////////////////////////////////////////////////////////////////////////////////////
    // In SMASH, one event can have many ensembles. Each of these ensembles is a separate
    // entry in the ROOT file. Moreover, one event can have multiple output times, leading
    // to even more entries. This needs to be appropriately handled. 

    /////////////////////////////////////////
    // If only_flow_at_final_output == true, only calculate flow at simulation end time
    //std::cout << "current_t = " << ROOT_file->current_t << std::endl;
    if (only_flow_at_final_output_ == true) {    
      if ( ROOT_file->current_t < config_info.End_Time() ) {
	continue;
      }
    }
    
    /////////////////////////////////////////
    // Establish which time step within an event this entry corresponds to   
    const int time_step_index = only_flow_at_final_output_ ?
      0 : (i_tree_entry % entries_per_event ) / ROOT_file->n_ensembles();

    /////////////////////////////////////////
    // Get info on time from the tree
    const double time = ROOT_file->t[0];

    
    
    //////////////////////////////////////////////////////////////////////////////////////
    // A loop over all particles in the entry, with if statements in the order from the
    // highest multiplicity to the lowest (approximately, and based on expectations for
    // FXT energies) 
    
    for (int index = 0; index < ROOT_file->npart; index++) {
      // Use switch logic to speed up
      switch ( ROOT_file->pdgcode[index] ) {
	///////////////////////////////////////
	// protons
      case 2212: {
	// count protons
	//proton_count[time_step_index]++;
	// add to N vs. y histograms
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_proton_in_4pi_[time_step_index],
		    proton_pT_min_, proton_pT_max_);
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_proton_at_mid_y_[time_step_index],
		    proton_pT_min_, proton_pT_max_);
	// add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_proton_basic_[time_step_index],
		      p_v2_proton_basic_[time_step_index],
		      p_v3_proton_basic_[time_step_index],
		      proton_pT_min_, proton_pT_max_);
	fill_integrated_v1_v2_v3
	  (ROOT_file->p0[index],
	   ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
	   p_v1_proton_basic_integrated_in_4pi_[time_step_index],
	   p_v2_proton_basic_integrated_in_4pi_[time_step_index],
	   p_v3_proton_basic_integrated_in_4pi_[time_step_index],
	   proton_pT_min_, proton_pT_max_);
	fill_integrated_v1_v2_v3
	  (ROOT_file->p0[index],
	   ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
	   p_v1_proton_basic_integrated_at_mid_y_[time_step_index],
	   p_v2_proton_basic_integrated_at_mid_y_[time_step_index],
	   p_v3_proton_basic_integrated_at_mid_y_[time_step_index],
	   proton_pT_min_, proton_pT_max_);
	
	break;
      }
	///////////////////////////////////////
	// pi minuses
      case -211: {
	// count pi minuses
	pi_minus_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_pi_minus_in_4pi_[time_step_index], pion_pT_min_, pion_pT_max_);
	// add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_pi_minus_basic_[time_step_index],
		      p_v2_pi_minus_basic_[time_step_index],
		      p_v3_pi_minus_basic_[time_step_index],
		      pion_pT_min_, pion_pT_max_);
	break;
      }
	///////////////////////////////////////
	// pi pluses
      case 211: {
	// count pi pluses
	pi_plus_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_pi_plus_in_4pi_[time_step_index], pion_pT_min_, pion_pT_max_);
	// add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_pi_plus_basic_[time_step_index],
		      p_v2_pi_plus_basic_[time_step_index],
		      p_v3_pi_plus_basic_[time_step_index],
		      pion_pT_min_, pion_pT_max_);
	break;
      }
	///////////////////////////////////////
	// deuterons
      case 1000010020: {
	// count deuterons
        deuteron_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_deuteron_in_4pi_[time_step_index],
		    deuteron_pT_min_, deuteron_pT_max_);
        // add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_deuteron_basic_[time_step_index],
		      p_v2_deuteron_basic_[time_step_index],
		      p_v3_deuteron_basic_[time_step_index],
		      deuteron_pT_min_, deuteron_pT_max_);
	break;
      }
	///////////////////////////////////////
	// K pluses
      case 321: {
	// count K pluses
	kaon_plus_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_kaon_plus_in_4pi_[time_step_index],
		    kaon_pT_min_, kaon_pT_max_);
	// add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_kaon_plus_basic_[time_step_index],
		      p_v2_kaon_plus_basic_[time_step_index],
		      p_v3_kaon_plus_basic_[time_step_index],
		      kaon_pT_min_, kaon_pT_max_);
	break;
      }
	///////////////////////////////////////
	// lambdas
      case 3122: {
	// count lambdas
        lambda_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_lambda_in_4pi_[time_step_index],
		    lambda_pT_min_, lambda_pT_max_);
        // add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_lambda_basic_[time_step_index],
		      p_v2_lambda_basic_[time_step_index],
		      p_v3_lambda_basic_[time_step_index],
		      lambda_pT_min_, lambda_pT_max_);
	break;
      }
	///////////////////////////////////////
	// K minuses
      case -321: {
	// count K minuses
	kaon_minus_count[time_step_index]++;
	// add to dN/dy histogram
	fill_N_vs_y(ROOT_file->p0[index],
		    ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		    h_N_kaon_minus_in_4pi_[time_step_index],
		    kaon_pT_min_, kaon_pT_max_);
	// add to flow TProfiles
	fill_v1_v2_v3(ROOT_file->p0[index],
		      ROOT_file->px[index], ROOT_file->py[index], ROOT_file->pz[index],
		      p_v1_kaon_minus_basic_[time_step_index],
		      p_v2_kaon_minus_basic_[time_step_index],
		      p_v3_kaon_minus_basic_[time_step_index],
		      kaon_pT_min_, kaon_pT_max_);
	break;
      }
	
      }

    } // for (int index = 0; index < npart; index++) {

  } // for (int i_tree_entry = 0; i_tree_entry < n_entries; i_tree_entry++) {

  std::cout << "\n\nFinished filling basic flow histograms \n" << std::endl;



  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Fit v1(y), v2(y), and v3(y) of protons:
  // p_v1_proton_basic, p_v2_proton_basic, p_v3_proton_basic
  // and pions:
  // p_v1_pi_plus_basic, p_v2_pi_plus_basic, p_v3_pi_plus_basic,
  // p_v1_pi_minus_basic, p_v2_pi_minus_basic, p_v3_pi_minus_basic
  //
  // Also, save data
  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Fitting and saving basic flow results               :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n"
	    << color::RESET << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Fit and plot flow results
  for (int i = 0; i < number_of_time_steps_; i++) {
    // Fit flow fata only for the selected times, if provided; otherwise, fit flow data
    // for all possible output times
    if ( !output_times.empty() && !output_times.count(ROOT_file->current_t_steps()[i]) ) {
      continue;
    }

    // Check that the integrated vn TProfiles work:
    double total_entries_in_4pi{};
    double v1_proton_in_4pi{};
    double v2_proton_in_4pi{};
    double v3_proton_in_4pi{};
    double total_entries_at_mid_y{};
    double v1_proton_at_mid_y{};
    double v2_proton_at_mid_y{};
    double v3_proton_at_mid_y{};
    // Loop over all rapidity bins
    for (int j = 1; j < (y_number_of_bins_ + 1); j++) {
      // Access N histogram at j-th bin
      const double h_N_proton_in_4pi_bin_content =
	h_N_proton_in_4pi_[i]->GetBinContent(j);
      const double sgn_y = (p_v1_proton_basic_[i]->GetBinCenter(j) > 0) ? 1.0 : -1.0;
      // Acces vn TProfiles at j-th bin
      const double p_v1_proton_bin_content = p_v1_proton_basic_[i]->GetBinContent(j);
      const double p_v2_proton_bin_content = p_v2_proton_basic_[i]->GetBinContent(j);
      const double p_v3_proton_bin_content = p_v3_proton_basic_[i]->GetBinContent(j);
      // Add to 4 pi sums
      total_entries_in_4pi += h_N_proton_in_4pi_bin_content;
      v1_proton_in_4pi += sgn_y * h_N_proton_in_4pi_bin_content * p_v1_proton_bin_content;
      v2_proton_in_4pi += h_N_proton_in_4pi_bin_content * p_v2_proton_bin_content;
      v3_proton_in_4pi += sgn_y * h_N_proton_in_4pi_bin_content * p_v3_proton_bin_content;
      // Add to midrapidity sums, if applicable:
      // Get bin witdh and bin center
      const double h_N_proton_bin_width = h_N_proton_in_4pi_[i]->GetBinWidth(j);
      const double h_N_proton_bin_center = h_N_proton_in_4pi_[i]->GetBinCenter(j);
      const double lower_edge = h_N_proton_bin_center - 0.5 * h_N_proton_bin_width;
      const double upper_edge = h_N_proton_bin_center + 0.5 * h_N_proton_bin_width;
      if ( (lower_edge >= y_mid_min_) && (upper_edge <= y_mid_max_) ) {
	total_entries_at_mid_y += h_N_proton_in_4pi_bin_content;
	v1_proton_at_mid_y +=
	  sgn_y * h_N_proton_in_4pi_bin_content * p_v1_proton_bin_content;
	v2_proton_at_mid_y += h_N_proton_in_4pi_bin_content * p_v2_proton_bin_content;
	v3_proton_at_mid_y +=
	  sgn_y * h_N_proton_in_4pi_bin_content * p_v3_proton_bin_content;
      }
    }
    // Normalize
    v1_proton_in_4pi /= total_entries_in_4pi;
    v2_proton_in_4pi /= total_entries_in_4pi;
    v3_proton_in_4pi /= total_entries_in_4pi;
    v1_proton_at_mid_y /= total_entries_at_mid_y;
    v2_proton_at_mid_y /= total_entries_at_mid_y;
    v3_proton_at_mid_y /= total_entries_at_mid_y;

    if (cfg.verbose) {
      std::cout << "\n                          v1_proton_in_4pi = " << v1_proton_in_4pi
		<< "\n                          v2_proton_in_4pi = " << v2_proton_in_4pi
		<< "\n                          v3_proton_in_4pi = " << v3_proton_in_4pi
		<< "\n"
		<< "\n   p_v1_proton_basic_integrated_in_4pi_[i] = "
		<< p_v1_proton_basic_integrated_in_4pi_[i]->GetBinContent(1)
		<< "\n   p_v2_proton_basic_integrated_in_4pi_[i] = "
		<< p_v2_proton_basic_integrated_in_4pi_[i]->GetBinContent(1)
		<< "\n   p_v3_proton_basic_integrated_in_4pi_[i] = "
		<< p_v3_proton_basic_integrated_in_4pi_[i]->GetBinContent(1)
		<< "\n"
		<< "\n                        v1_proton_at_mid_y = " << v1_proton_at_mid_y
		<< "\n                        v2_proton_at_mid_y = " << v2_proton_at_mid_y
		<< "\n                        v3_proton_at_mid_y = " << v3_proton_at_mid_y
		<< "\n"
		<< "\n p_v1_proton_basic_integrated_at_mid_y_[i] = "
		<< p_v1_proton_basic_integrated_at_mid_y_[i]->GetBinContent(1)
		<< "\n p_v2_proton_basic_integrated_at_mid_y_[i] = "
		<< p_v2_proton_basic_integrated_at_mid_y_[i]->GetBinContent(1)
		<< "\n p_v3_proton_basic_integrated_at_mid_y_[i] = "
		<< p_v3_proton_basic_integrated_at_mid_y_[i]->GetBinContent(1)
		<< std::endl;

    
      std::cout << "\n # of v1 bins at mid y = "
		<< p_v1_proton_basic_integrated_at_mid_y_[i]->GetNbinsX()
		<< "\n # of v2 bins at mid y = "
		<< p_v2_proton_basic_integrated_at_mid_y_[i]->GetNbinsX()
		<< "\n # of v3 bins at mid y = "
		<< p_v3_proton_basic_integrated_at_mid_y_[i]->GetNbinsX()
		<< std::endl;
    }

    // protons:
    fit_v1_v2_v3_make_plots_save_fit_data
      (p_v1_proton_basic_[i], p_v2_proton_basic_[i], p_v3_proton_basic_[i],
       ROOT_file->n_events(), config_info.Range_or_Value(), ROOT_file->time_steps()[i],
       cfg, true);
    fit_v1_v2_v3_make_plots_save_fit_data
      (p_v1_proton_basic_integrated_in_4pi_[i],
       p_v2_proton_basic_integrated_in_4pi_[i],
       p_v3_proton_basic_integrated_in_4pi_[i],
       ROOT_file->n_events(), config_info.Range_or_Value(), ROOT_file->time_steps()[i],
       cfg, true);
    fit_v1_v2_v3_make_plots_save_fit_data
      (p_v1_proton_basic_integrated_at_mid_y_[i],
       p_v2_proton_basic_integrated_at_mid_y_[i],
       p_v3_proton_basic_integrated_at_mid_y_[i],
       ROOT_file->n_events(), config_info.Range_or_Value(), ROOT_file->time_steps()[i],
       cfg, true);
    // pi pluses:
    fit_v1_v2_v3_make_plots_save_fit_data
      (p_v1_pi_plus_basic_[i], p_v2_pi_plus_basic_[i], p_v3_pi_plus_basic_[i],
       ROOT_file->n_events(), config_info.Range_or_Value(), ROOT_file->time_steps()[i],
       cfg, true);
    // pi minuses:
    fit_v1_v2_v3_make_plots_save_fit_data
      (p_v1_pi_minus_basic_[i], p_v2_pi_minus_basic_[i], p_v3_pi_minus_basic_[i],
       ROOT_file->n_events(), config_info.Range_or_Value(), ROOT_file->time_steps()[i],
       cfg, true);
  }

  std::cout << "Finished fitting flow histograms" << std::endl;
  

  
  ///////////////////////////////////////////
  // Save flow vs. rapidty data    
  for (int i = 0; i < number_of_time_steps_; i++) {
    // Save data only for the selected times, if provided; otherwise, save flow data for
    // all possible output times
    if ( !output_times.empty() && !output_times.count(ROOT_file->current_t_steps()[i]) ) {
      continue;
    }

    save_flow_vs_rapidity_data
      (p_v1_proton_basic_[i], p_v1_deuteron_basic_[i], p_v1_lambda_basic_[i],
       p_v1_pi_plus_basic_[i], p_v1_pi_minus_basic_[i], p_v1_kaon_plus_basic_[i],
       p_v1_kaon_minus_basic_[i],
       p_v2_proton_basic_[i], p_v2_deuteron_basic_[i], p_v2_lambda_basic_[i],
       p_v2_pi_plus_basic_[i], p_v2_pi_minus_basic_[i], p_v2_kaon_plus_basic_[i],
       p_v2_kaon_minus_basic_[i],
       p_v3_proton_basic_[i], p_v3_deuteron_basic_[i], p_v3_lambda_basic_[i],
       p_v3_pi_plus_basic_[i], p_v3_pi_minus_basic_[i], p_v3_kaon_plus_basic_[i],
       p_v3_kaon_minus_basic_[i],
       ROOT_file->n_events(),
       config_info.Range_or_Value(),
       ROOT_file->time_steps()[i]);
  }   

  

  ///////////////////////////////////////////
  // Print and save yield data
  const int real_event_equivalent =
    ROOT_file->n_events() * ROOT_file->n_test() * ROOT_file->n_ensembles();

  save_yield_data(proton_count, deuteron_count, lambda_count,
		  pi_plus_count, pi_minus_count, kaon_plus_count, kaon_minus_count,   
		  ROOT_file->current_t_steps(), output_times,
		  p_v1_proton_basic_[0], ROOT_file->n_events(), real_event_equivalent,
		  cfg);
  

  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Update the control variable and exit
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  // Re-enable all branches
  ROOT_file->fChain->SetBranchStatus("*", 1); 
    
  basic_flow_analysis_done_ = true;

  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Finished basic flow analysis                        :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n\n"
	    << color::RESET << std::endl;
  
}



void Flow::basic_flow_time_evolution_binned_in_y
  (const std::unique_ptr<ReadParticles>& ROOT_file, Config cfg) {

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform basic flow time evolution analysis
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Starting basic flow time evolution analysis in y bins :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ."
	    << color::RESET << std::endl;

  ///////////////////////////////////////////
  // Basic checks
  if ( only_flow_at_final_output_ ) {
    throw std::runtime_error("Need to initialize ColliderFlow object with multiple time "
			     "steps to study the time evolution");
  }
  if ( !basic_flow_analysis_done_ ) {
    throw std::runtime_error("Need to perform basic_flow() analysis first to study "
			     "the time evolution of basic flow");
  }


  
  ///////////////////////////////////////////
  // Define variables:
  // time step
  const double dt = ROOT_file->output_interval();
  // previous vn
  std::vector<double> v1_proton_basic_previous_value(y_number_of_bins_, 0.0);
  std::vector<double> v2_proton_basic_previous_value(y_number_of_bins_, 0.0);
  std::vector<double> v3_proton_basic_previous_value(y_number_of_bins_, 0.0);
  std::vector<double> v1_proton_basic_previous_error(y_number_of_bins_, 0.0);
  std::vector<double> v2_proton_basic_previous_error(y_number_of_bins_, 0.0);
  std::vector<double> v3_proton_basic_previous_error(y_number_of_bins_, 0.0);
  // current vn
  std::vector<double> v1_proton_basic_current_value(y_number_of_bins_, 0.0);
  std::vector<double> v2_proton_basic_current_value(y_number_of_bins_, 0.0);
  std::vector<double> v3_proton_basic_current_value(y_number_of_bins_, 0.0);
  std::vector<double> v1_proton_basic_current_error(y_number_of_bins_, 0.0);
  std::vector<double> v2_proton_basic_current_error(y_number_of_bins_, 0.0);
  std::vector<double> v3_proton_basic_current_error(y_number_of_bins_, 0.0);

  // Variables to keep track of MF and collision contributions through values of current_t
  double current_t{};
  double previous_current_t{};

  // Auxiliary variables to store the d_vn{coll,MF,total} contributions (differential
  // changes in vn due to {coll,MF,total}) for a given current_t, at all bins:  
  std::vector<double> dv1coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv2coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv3coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv1coll_proton_basic_err(y_number_of_bins_, 0.0);
  std::vector<double> dv2coll_proton_basic_err(y_number_of_bins_, 0.0);
  std::vector<double> dv3coll_proton_basic_err(y_number_of_bins_, 0.0);
  std::vector<double> dv1MF_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv2MF_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv3MF_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> dv1MF_proton_basic_err(y_number_of_bins_, 0.0);
  std::vector<double> dv2MF_proton_basic_err(y_number_of_bins_, 0.0);
  std::vector<double> dv3MF_proton_basic_err(y_number_of_bins_, 0.0);
  // Auxiliary variables to store integrated contributions to vn{coll,MF,total}:
  std::vector<double> int_dv1coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> int_dv2coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> int_dv3coll_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> int_dv1MF_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> int_dv2MF_proton_basic(y_number_of_bins_, 0.0);
  std::vector<double> int_dv3MF_proton_basic(y_number_of_bins_, 0.0);

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Loop over TProfiles and extract differential increases in vnMF, vncoll, and vntotal,
  // then construct time derivatives, then construct vnMF, vncoll, and vntotal

  // We keep track of the index of TGraph entries with these variables; they start from -1
  // because they are updated at the top of the loops (the first update then immediately
  // gives i_coll_entry = 0, i_MF_entry = 0)
  int i_coll_entry = -1;
  int i_MF_entry = -1;
  
  for (int i = 0; i < number_of_time_steps_; i++) {
    /////////////////////////////////////////
    // Get previous and current values of current_t:
    previous_current_t = current_t;
    current_t = ROOT_file->current_t_steps()[i];
    
    for (int j = 0; j < y_number_of_bins_; j++) {
      ///////////////////////////////////////
      // Get previous and current values of v_n(y):      
      // Previous values:
      v1_proton_basic_previous_value[j] = v1_proton_basic_current_value[j];
      v1_proton_basic_previous_error[j] = v1_proton_basic_current_error[j];
      v2_proton_basic_previous_value[j] = v2_proton_basic_current_value[j];
      v2_proton_basic_previous_error[j] = v2_proton_basic_current_error[j];
      v3_proton_basic_previous_value[j] = v3_proton_basic_current_value[j];      
      v3_proton_basic_previous_error[j] = v3_proton_basic_current_error[j];
      // Current values:
      // Histograms and TProfiles start bins from 1, not from 0
      const int j_bin = j + 1;
      v1_proton_basic_current_value[j] = p_v1_proton_basic_[i]->GetBinContent(j_bin);     
      v1_proton_basic_current_error[j] = p_v1_proton_basic_[i]->GetBinError(j_bin);
      v2_proton_basic_current_value[j] = p_v2_proton_basic_[i]->GetBinContent(j_bin);     
      v2_proton_basic_current_error[j] = p_v2_proton_basic_[i]->GetBinError(j_bin);
      v3_proton_basic_current_value[j] = p_v3_proton_basic_[i]->GetBinContent(j_bin);     
      v3_proton_basic_current_error[j] = p_v3_proton_basic_[i]->GetBinError(j_bin);

      ///////////////////////////////////////
      // Compute time differentials and assign derivatives to a correct TGraph
      if ( i == 0 ) {
	// Cannot compute the time differential as there is no previous step
	continue;
      }


      
      ///////////////////////////////////////
      // We know the structure of the file: if current_t > previous_current_t, the output
      // is due to the collision contribution, while for current_t = previous_current_t
      // the output is due to the MF contribution; at current_t = previous_current_t, one
      // can also record the total contribution.

      if (current_t > previous_current_t) {
	// If j == 0, update the TGraph entry index
	if ( j == 0) {
	  i_coll_entry++;
	}
	
	/////////////////////////////////////
	// Compute observables:
	// vncoll
	dv1coll_proton_basic[j] =
	  v1_proton_basic_current_value[j] - v1_proton_basic_previous_value[j];
	dv2coll_proton_basic[j] =
	  v2_proton_basic_current_value[j] - v2_proton_basic_previous_value[j];
	dv3coll_proton_basic[j] =
	  v3_proton_basic_current_value[j] - v3_proton_basic_previous_value[j];
	// vncoll error
	dv1coll_proton_basic_err[j] =
	  std::sqrt( std::pow(v1_proton_basic_current_error[j], 2.0)
		     + std::pow(v1_proton_basic_previous_error[j], 2.0) );
	dv2coll_proton_basic_err[j] =
	  std::sqrt( std::pow(v2_proton_basic_current_error[j], 2.0)
		     + std::pow(v2_proton_basic_previous_error[j], 2.0) );
	dv3coll_proton_basic_err[j] =
	  std::sqrt( std::pow(v3_proton_basic_current_error[j], 2.0)
		     + std::pow(v3_proton_basic_previous_error[j], 2.0) );
	// Integrate to obtain int dvncoll
	int_dv1coll_proton_basic[j] += dv1coll_proton_basic[j];
	int_dv2coll_proton_basic[j] += dv2coll_proton_basic[j];
	int_dv3coll_proton_basic[j] += dv3coll_proton_basic[j];

	/////////////////////////////////////
	// Record dvncoll and int dvncoll
	// dv1coll/dt
	g_dv1coll_dt_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, dv1coll_proton_basic[j]/dt);
	g_dv1coll_dt_proton_basic_[j]->
	  SetPointError(i_coll_entry, 0.0, dv1coll_proton_basic_err[j]/dt);
	// dv2coll/dt
	g_dv2coll_dt_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, dv2coll_proton_basic[j]/dt);
	g_dv2coll_dt_proton_basic_[j]->
	  SetPointError(i_coll_entry, 0.0, dv2coll_proton_basic_err[j]/dt);
	// dv3coll/dt
	g_dv3coll_dt_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, dv3coll_proton_basic[j]/dt);
	g_dv3coll_dt_proton_basic_[j]->SetPointError
	  (i_coll_entry, 0.0, dv3coll_proton_basic_err[j]/dt);
	// int dvncoll
	g_v1coll_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, int_dv1coll_proton_basic[j]);	
	g_v2coll_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, int_dv2coll_proton_basic[j]);	
	g_v3coll_proton_basic_[j]->
	  SetPoint(i_coll_entry, current_t, int_dv3coll_proton_basic[j]);
      } else if (current_t == previous_current_t) {
	// If j == 0, update the TGraph entry index
	if ( j == 0) {
	  i_MF_entry++;
	}

	/////////////////////////////////////
	// Compute observables:
	// vnMF
	dv1MF_proton_basic[j] =
	  v1_proton_basic_current_value[j] - v1_proton_basic_previous_value[j];
	dv2MF_proton_basic[j] =
	  v2_proton_basic_current_value[j] - v2_proton_basic_previous_value[j];
	dv3MF_proton_basic[j] =
	  v3_proton_basic_current_value[j] - v3_proton_basic_previous_value[j];
	// vnMF error
	dv1MF_proton_basic_err[j] =
	  std::sqrt( std::pow(v1_proton_basic_current_error[j], 2.0)
		     + std::pow(v1_proton_basic_previous_error[j], 2.0) );
	dv2MF_proton_basic_err[j] =
	  std::sqrt( std::pow(v2_proton_basic_current_error[j], 2.0)
		     + std::pow(v2_proton_basic_previous_error[j], 2.0) );
	dv3MF_proton_basic_err[j] =
	  std::sqrt( std::pow(v3_proton_basic_current_error[j], 2.0)
		     + std::pow(v3_proton_basic_previous_error[j], 2.0) );
	// dv1total/dt
	const double dv1total = dv1coll_proton_basic[j] + dv1MF_proton_basic[j];
	const double dv1total_err =
	  std::sqrt( std::pow(dv1coll_proton_basic_err[j], 2.0) +
		     std::pow(dv1MF_proton_basic_err[j], 2.0) );
	// dv2total/dt
        const double dv2total = dv2coll_proton_basic[j] + dv2MF_proton_basic[j];
	const double dv2total_err =
	  std::sqrt( std::pow(dv2coll_proton_basic_err[j], 2.0) +
		     std::pow(dv2MF_proton_basic_err[j], 2.0) );
	// dv3total/dt
        const double dv3total = dv3coll_proton_basic[j] + dv3MF_proton_basic[j];
	const double dv3total_err =
	  std::sqrt( std::pow(dv3coll_proton_basic_err[j], 2.0) +
		     std::pow(dv3MF_proton_basic_err[j], 2.0) );
	// Integrate to obtain int dvnMF (and enable int dvntotal)
        int_dv1MF_proton_basic[j] += dv1MF_proton_basic[j];
	int_dv2MF_proton_basic[j] += dv2MF_proton_basic[j];
	int_dv3MF_proton_basic[j] += dv3MF_proton_basic[j];

	/////////////////////////////////////
	// Record dvnMF and int dvnMF, and dvntotal and int dvntotal:
	// dv1MF/dt
	g_dv1MF_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv1MF_proton_basic[j]/dt);
	g_dv1MF_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv1MF_proton_basic_err[j]/dt);
	// dv2MF/dt
	g_dv2MF_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv2MF_proton_basic[j]/dt);
	g_dv2MF_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv2MF_proton_basic_err[j]/dt);
	// dv3MF/dt
	g_dv3MF_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv3MF_proton_basic[j]/dt);
	g_dv3MF_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv3MF_proton_basic_err[j]/dt);
	// dv1total/dt
        g_dv1total_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv1total/dt);
	g_dv1total_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv1total_err/dt);
	// dv2total/dt
	g_dv2total_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv2total/dt);
	g_dv2total_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv2total_err/dt);
	// dv3total/dt
	g_dv3total_dt_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, dv3total/dt);
	g_dv3total_dt_proton_basic_[j]->
	  SetPointError(i_MF_entry, 0.0, dv3total_err/dt);
	// int dvnMF
	g_v1MF_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, int_dv1MF_proton_basic[j]);	
	g_v2MF_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, int_dv2MF_proton_basic[j]);	
	g_v3MF_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t, int_dv3MF_proton_basic[j]);
	// int dvntotal
	g_v1total_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t,
		   int_dv1coll_proton_basic[j] + int_dv1MF_proton_basic[j]);
	g_v2total_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t,
		   int_dv2coll_proton_basic[j] + int_dv2MF_proton_basic[j]);
	g_v3total_proton_basic_[j]->
	  SetPoint(i_MF_entry, current_t,
		   int_dv3coll_proton_basic[j] + int_dv3MF_proton_basic[j]);
      }

    }
    
  }



  ////////////////////////////////////////////////////////////////////////////////////////
  // Plot and save data on dvn/dt

  const int real_event_equivalent =
    ROOT_file->n_events() * ROOT_file->n_test() * ROOT_file->n_ensembles();

  const int every_nth_entry = 5;
  
  // Loop over all bins and plot observables
  for (int i = 0; i < y_number_of_bins_; i++) {
    plot_and_save_flow_time_evolution_data
      (g_dv1coll_dt_proton_basic_[i],
       g_dv1MF_dt_proton_basic_[i],
       g_dv1total_dt_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.15, 0.15,
       every_nth_entry, cfg, true);

    plot_and_save_flow_time_evolution_data
      (g_dv2coll_dt_proton_basic_[i],
       g_dv2MF_dt_proton_basic_[i],
       g_dv2total_dt_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.05, 0.05,
       every_nth_entry, cfg, true);

    plot_and_save_flow_time_evolution_data
      (g_dv3coll_dt_proton_basic_[i],
       g_dv3MF_dt_proton_basic_[i],
       g_dv3total_dt_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.075, 0.075,
       every_nth_entry, cfg, true);

    plot_and_save_flow_time_evolution_data
      (g_v1coll_proton_basic_[i],
       g_v1MF_proton_basic_[i],
       g_v1total_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.15, 0.15,
       every_nth_entry, cfg, true);

    plot_and_save_flow_time_evolution_data
      (g_v2coll_proton_basic_[i],
       g_v2MF_proton_basic_[i],
       g_v2total_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.1, 0.1,
       every_nth_entry, cfg, true);

    plot_and_save_flow_time_evolution_data
      (g_v3coll_proton_basic_[i],
       g_v3MF_proton_basic_[i],
       g_v3total_proton_basic_[i],
       ROOT_file->n_events(), real_event_equivalent,
       ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
       // these values are just guesses
       -0.075, 0.075,
       every_nth_entry, cfg, true);
  }

  

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Finished basic flow time evolution analysis in y bins  :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n\n"
	    << color::RESET << std::endl;
    
}




void Flow::basic_flow_time_evolution_in_4pi_and_at_mid_y
  (const std::unique_ptr<ReadParticles>& ROOT_file, Config cfg) {

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform basic flow time evolution analysis
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Starting basic flow time evolution analysis in 4pi and at mid y :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << color::RESET << std::endl;

  ///////////////////////////////////////////
  // Basic checks
  if ( only_flow_at_final_output_ ) {
    throw std::runtime_error("Need to initialize ColliderFlow object with multiple time "
			     "steps to study the time evolution");
  }
  if ( !basic_flow_analysis_done_ ) {
    throw std::runtime_error("Need to perform basic_flow() analysis first to study "
			     "the time evolution of basic flow");
  }


  
  ///////////////////////////////////////////
  // Define variables:
  // time step
  const double dt = ROOT_file->output_interval();
  // previous vn
  double v1_proton_basic_in_4pi_previous_value{};
  double v2_proton_basic_in_4pi_previous_value{};
  double v3_proton_basic_in_4pi_previous_value{};
  double v1_proton_basic_in_4pi_previous_error{};
  double v2_proton_basic_in_4pi_previous_error{};
  double v3_proton_basic_in_4pi_previous_error{};
  double v1_proton_basic_at_mid_y_previous_value{};
  double v2_proton_basic_at_mid_y_previous_value{};
  double v3_proton_basic_at_mid_y_previous_value{};
  double v1_proton_basic_at_mid_y_previous_error{};
  double v2_proton_basic_at_mid_y_previous_error{};
  double v3_proton_basic_at_mid_y_previous_error{};
  // current vn
  double v1_proton_basic_in_4pi_current_value{};
  double v2_proton_basic_in_4pi_current_value{};
  double v3_proton_basic_in_4pi_current_value{};
  double v1_proton_basic_in_4pi_current_error{};
  double v2_proton_basic_in_4pi_current_error{};
  double v3_proton_basic_in_4pi_current_error{};
  double v1_proton_basic_at_mid_y_current_value{};
  double v2_proton_basic_at_mid_y_current_value{};
  double v3_proton_basic_at_mid_y_current_value{};
  double v1_proton_basic_at_mid_y_current_error{};
  double v2_proton_basic_at_mid_y_current_error{};
  double v3_proton_basic_at_mid_y_current_error{};

  // Variables to keep track of MF and collision contributions through values of current_t
  double current_t{};
  double previous_current_t{};

  // Auxiliary variables to store the d_vn{coll,MF,total} contributions (differential
  // changes in vn due to {coll,MF,total}) for a given current_t:  
  double dv1coll_proton_basic_in_4pi{};
  double dv2coll_proton_basic_in_4pi{};
  double dv3coll_proton_basic_in_4pi{};
  double dv1coll_proton_basic_in_4pi_err{};
  double dv2coll_proton_basic_in_4pi_err{};
  double dv3coll_proton_basic_in_4pi_err{};
  double dv1MF_proton_basic_in_4pi{};
  double dv2MF_proton_basic_in_4pi{};
  double dv3MF_proton_basic_in_4pi{};
  double dv1MF_proton_basic_in_4pi_err{};
  double dv2MF_proton_basic_in_4pi_err{};
  double dv3MF_proton_basic_in_4pi_err{};
  double dv1coll_proton_basic_at_mid_y{};
  double dv2coll_proton_basic_at_mid_y{};
  double dv3coll_proton_basic_at_mid_y{};
  double dv1coll_proton_basic_at_mid_y_err{};
  double dv2coll_proton_basic_at_mid_y_err{};
  double dv3coll_proton_basic_at_mid_y_err{};
  double dv1MF_proton_basic_at_mid_y{};
  double dv2MF_proton_basic_at_mid_y{};
  double dv3MF_proton_basic_at_mid_y{};
  double dv1MF_proton_basic_at_mid_y_err{};
  double dv2MF_proton_basic_at_mid_y_err{};
  double dv3MF_proton_basic_at_mid_y_err{};
  // Auxiliary variables to store integrated contributions to vn{coll,MF,total}:
  double int_dv1coll_proton_basic_in_4pi{};
  double int_dv2coll_proton_basic_in_4pi{};
  double int_dv3coll_proton_basic_in_4pi{};
  double int_dv1MF_proton_basic_in_4pi{};
  double int_dv2MF_proton_basic_in_4pi{};
  double int_dv3MF_proton_basic_in_4pi{};
  double int_dv1coll_proton_basic_at_mid_y{};
  double int_dv2coll_proton_basic_at_mid_y{};
  double int_dv3coll_proton_basic_at_mid_y{};
  double int_dv1MF_proton_basic_at_mid_y{};
  double int_dv2MF_proton_basic_at_mid_y{};
  double int_dv3MF_proton_basic_at_mid_y{};

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Loop over TProfiles and extract differential increases in vnMF, vncoll, and vntotal,
  // then construct time derivatives, then construct vnMF, vncoll, and vntotal

  // We keep track of the index of TGraph entries with these variables; they start from -1
  // because they are updated at the top of the loops (the first update then immediately
  // gives i_coll_entry = 0, i_MF_entry = 0)
  int i_coll_entry = -1;
  int i_MF_entry = -1;
  
  for (int i = 0; i < number_of_time_steps_; i++) {
    /////////////////////////////////////////
    // Get previous and current values of current_t:
    previous_current_t = current_t;
    current_t = ROOT_file->current_t_steps()[i];
    

    ///////////////////////////////////////
    // Get previous and current values of v_n in 4pi and at mid y:      
    // Previous values:
    // in 4pi
    v1_proton_basic_in_4pi_previous_value = v1_proton_basic_in_4pi_current_value;
    v1_proton_basic_in_4pi_previous_error = v1_proton_basic_in_4pi_current_error;
    v2_proton_basic_in_4pi_previous_value = v2_proton_basic_in_4pi_current_value;
    v2_proton_basic_in_4pi_previous_error = v2_proton_basic_in_4pi_current_error;
    v3_proton_basic_in_4pi_previous_value = v3_proton_basic_in_4pi_current_value;      
    v3_proton_basic_in_4pi_previous_error = v3_proton_basic_in_4pi_current_error;
    // at mid y
    v1_proton_basic_at_mid_y_previous_value = v1_proton_basic_at_mid_y_current_value;
    v1_proton_basic_at_mid_y_previous_error = v1_proton_basic_at_mid_y_current_error;
    v2_proton_basic_at_mid_y_previous_value = v2_proton_basic_at_mid_y_current_value;
    v2_proton_basic_at_mid_y_previous_error = v2_proton_basic_at_mid_y_current_error;
    v3_proton_basic_at_mid_y_previous_value = v3_proton_basic_at_mid_y_current_value;      
    v3_proton_basic_at_mid_y_previous_error = v3_proton_basic_at_mid_y_current_error;
    
    // Current values:
    // Histograms and TProfiles start bins from 1, not from 0;
    // for 4pi and at mid y analysis, there is only one bin so we hard-code j_bin to 1
    const int j_bin = 1;
    // in 4pi
    v1_proton_basic_in_4pi_current_value =
      p_v1_proton_basic_integrated_in_4pi_[i]->GetBinContent(j_bin);     
    v1_proton_basic_in_4pi_current_error =
      p_v1_proton_basic_integrated_in_4pi_[i]->GetBinError(j_bin);
    v2_proton_basic_in_4pi_current_value =
      p_v2_proton_basic_integrated_in_4pi_[i]->GetBinContent(j_bin);     
    v2_proton_basic_in_4pi_current_error =
      p_v2_proton_basic_integrated_in_4pi_[i]->GetBinError(j_bin);
    v3_proton_basic_in_4pi_current_value =
      p_v3_proton_basic_integrated_in_4pi_[i]->GetBinContent(j_bin);     
    v3_proton_basic_in_4pi_current_error =
      p_v3_proton_basic_integrated_in_4pi_[i]->GetBinError(j_bin);
    // at mid y
    v1_proton_basic_at_mid_y_current_value =
      p_v1_proton_basic_integrated_at_mid_y_[i]->GetBinContent(j_bin);     
    v1_proton_basic_at_mid_y_current_error =
      p_v1_proton_basic_integrated_at_mid_y_[i]->GetBinError(j_bin);
    v2_proton_basic_at_mid_y_current_value =
      p_v2_proton_basic_integrated_at_mid_y_[i]->GetBinContent(j_bin);     
    v2_proton_basic_at_mid_y_current_error =
      p_v2_proton_basic_integrated_at_mid_y_[i]->GetBinError(j_bin);
    v3_proton_basic_at_mid_y_current_value =
      p_v3_proton_basic_integrated_at_mid_y_[i]->GetBinContent(j_bin);     
    v3_proton_basic_at_mid_y_current_error =
      p_v3_proton_basic_integrated_at_mid_y_[i]->GetBinError(j_bin);

    ///////////////////////////////////////
    // Compute time differentials and assign derivatives to a correct TGraph
    if ( i == 0 ) {
      // Cannot compute the time differential as there is no previous step
      continue;
    }



    ///////////////////////////////////////
    // We know the structure of the file: if current_t > previous_current_t, the output
    // is due to the collision contribution, while for current_t = previous_current_t
    // the output is due to the MF contribution; at current_t = previous_current_t, one
    // can also record the total contribution.

    if (current_t > previous_current_t) {
      // Update the TGraph entry index
      i_coll_entry++;
	
      /////////////////////////////////////
      // Compute observables:
      // vncoll in 4pi
      dv1coll_proton_basic_in_4pi =
	v1_proton_basic_in_4pi_current_value - v1_proton_basic_in_4pi_previous_value;
      dv2coll_proton_basic_in_4pi =
	v2_proton_basic_in_4pi_current_value - v2_proton_basic_in_4pi_previous_value;
      dv3coll_proton_basic_in_4pi =
	v3_proton_basic_in_4pi_current_value - v3_proton_basic_in_4pi_previous_value;
      // vncoll at mid y
      dv1coll_proton_basic_at_mid_y =
	v1_proton_basic_at_mid_y_current_value - v1_proton_basic_at_mid_y_previous_value;
      dv2coll_proton_basic_at_mid_y =
	v2_proton_basic_at_mid_y_current_value - v2_proton_basic_at_mid_y_previous_value;
      dv3coll_proton_basic_at_mid_y =
	v3_proton_basic_at_mid_y_current_value - v3_proton_basic_at_mid_y_previous_value;
      // vncoll error in 4 pi
      dv1coll_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v1_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v1_proton_basic_in_4pi_previous_error, 2.0) );
      dv2coll_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v2_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v2_proton_basic_in_4pi_previous_error, 2.0) );
      dv3coll_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v3_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v3_proton_basic_in_4pi_previous_error, 2.0) );
      // vncoll error at mid y
      dv1coll_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v1_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v1_proton_basic_at_mid_y_previous_error, 2.0) );
      dv2coll_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v2_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v2_proton_basic_at_mid_y_previous_error, 2.0) );
      dv3coll_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v3_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v3_proton_basic_at_mid_y_previous_error, 2.0) );
      // Integrate to obtain int dvncoll in 4pi
      int_dv1coll_proton_basic_in_4pi += dv1coll_proton_basic_in_4pi;
      int_dv2coll_proton_basic_in_4pi += dv2coll_proton_basic_in_4pi;
      int_dv3coll_proton_basic_in_4pi += dv3coll_proton_basic_in_4pi;
      // Integrate to obtain int dvncoll at mid y
      int_dv1coll_proton_basic_at_mid_y += dv1coll_proton_basic_at_mid_y;
      int_dv2coll_proton_basic_at_mid_y += dv2coll_proton_basic_at_mid_y;
      int_dv3coll_proton_basic_at_mid_y += dv3coll_proton_basic_at_mid_y;

      /////////////////////////////////////
      // Record dvncoll and int dvncoll
      // dv1coll/dt in 4pi
      g_dv1coll_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, dv1coll_proton_basic_in_4pi/dt);
      g_dv1coll_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_coll_entry, 0.0, dv1coll_proton_basic_in_4pi_err/dt);
      // dv1coll/dt at mid y
      g_dv1coll_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, dv1coll_proton_basic_at_mid_y/dt);
      g_dv1coll_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_coll_entry, 0.0, dv1coll_proton_basic_at_mid_y_err/dt);
      // dv2coll/dt in 4pi
      g_dv2coll_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, dv2coll_proton_basic_in_4pi/dt);
      g_dv2coll_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_coll_entry, 0.0, dv2coll_proton_basic_in_4pi_err/dt);
      // dv2coll/dt at mid y
      g_dv2coll_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, dv2coll_proton_basic_at_mid_y/dt);
      g_dv2coll_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_coll_entry, 0.0, dv2coll_proton_basic_at_mid_y_err/dt);
      // dv3coll/dt in 4pi
      g_dv3coll_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, dv3coll_proton_basic_in_4pi/dt);
      g_dv3coll_dt_proton_basic_integrated_in_4pi_->SetPointError
	(i_coll_entry, 0.0, dv3coll_proton_basic_in_4pi_err/dt);
      // dv3coll/dt at mid y
      g_dv3coll_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, dv3coll_proton_basic_at_mid_y/dt);
      g_dv3coll_dt_proton_basic_integrated_at_mid_y_->SetPointError
	(i_coll_entry, 0.0, dv3coll_proton_basic_at_mid_y_err/dt);
      // int dvncoll in 4pi
      g_v1coll_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, int_dv1coll_proton_basic_in_4pi);	
      g_v2coll_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, int_dv2coll_proton_basic_in_4pi);	
      g_v3coll_proton_basic_integrated_in_4pi_->
	SetPoint(i_coll_entry, current_t, int_dv3coll_proton_basic_in_4pi);
      // int dvncoll at mid y
      g_v1coll_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, int_dv1coll_proton_basic_at_mid_y);	
      g_v2coll_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, int_dv2coll_proton_basic_at_mid_y);	
      g_v3coll_proton_basic_integrated_at_mid_y_->
	SetPoint(i_coll_entry, current_t, int_dv3coll_proton_basic_at_mid_y);
      
    } else if (current_t == previous_current_t) {
      // Update the TGraph entry index
      i_MF_entry++;

      /////////////////////////////////////
      // Compute observables:
      // vnMF in 4pi
      dv1MF_proton_basic_in_4pi =
	v1_proton_basic_in_4pi_current_value - v1_proton_basic_in_4pi_previous_value;
      dv2MF_proton_basic_in_4pi =
	v2_proton_basic_in_4pi_current_value - v2_proton_basic_in_4pi_previous_value;
      dv3MF_proton_basic_in_4pi =
	v3_proton_basic_in_4pi_current_value - v3_proton_basic_in_4pi_previous_value;
      // vnMF at mid y
      dv1MF_proton_basic_at_mid_y =
	v1_proton_basic_at_mid_y_current_value - v1_proton_basic_at_mid_y_previous_value;
      dv2MF_proton_basic_at_mid_y =
	v2_proton_basic_at_mid_y_current_value - v2_proton_basic_at_mid_y_previous_value;
      dv3MF_proton_basic_at_mid_y =
	v3_proton_basic_at_mid_y_current_value - v3_proton_basic_at_mid_y_previous_value;
      // vnMF error in 4pi
      dv1MF_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v1_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v1_proton_basic_in_4pi_previous_error, 2.0) );
      dv2MF_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v2_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v2_proton_basic_in_4pi_previous_error, 2.0) );
      dv3MF_proton_basic_in_4pi_err =
	std::sqrt( std::pow(v3_proton_basic_in_4pi_current_error, 2.0)
		   + std::pow(v3_proton_basic_in_4pi_previous_error, 2.0) );
      // vnMF error at mid y
      dv1MF_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v1_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v1_proton_basic_at_mid_y_previous_error, 2.0) );
      dv2MF_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v2_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v2_proton_basic_at_mid_y_previous_error, 2.0) );
      dv3MF_proton_basic_at_mid_y_err =
	std::sqrt( std::pow(v3_proton_basic_at_mid_y_current_error, 2.0)
		   + std::pow(v3_proton_basic_at_mid_y_previous_error, 2.0) );
      // dv1total/dt in 4pi
      const double dv1total_in_4pi =
	dv1coll_proton_basic_in_4pi + dv1MF_proton_basic_in_4pi;
      const double dv1total_in_4pi_err =
	std::sqrt( std::pow(dv1coll_proton_basic_in_4pi_err, 2.0) +
		   std::pow(dv1MF_proton_basic_in_4pi_err, 2.0) );
      // dv1total/dt at mid y
      const double dv1total_at_mid_y =
	dv1coll_proton_basic_at_mid_y + dv1MF_proton_basic_at_mid_y;
      const double dv1total_at_mid_y_err =
	std::sqrt( std::pow(dv1coll_proton_basic_at_mid_y_err, 2.0) +
		   std::pow(dv1MF_proton_basic_at_mid_y_err, 2.0) );
      // dv2total/dt in 4pi
      const double dv2total_in_4pi =
	dv2coll_proton_basic_in_4pi + dv2MF_proton_basic_in_4pi;
      const double dv2total_in_4pi_err =
	std::sqrt( std::pow(dv2coll_proton_basic_in_4pi_err, 2.0) +
		   std::pow(dv2MF_proton_basic_in_4pi_err, 2.0) );
      // dv2total/dt at mid y
      const double dv2total_at_mid_y =
	dv2coll_proton_basic_at_mid_y + dv2MF_proton_basic_at_mid_y;
      const double dv2total_at_mid_y_err =
	std::sqrt( std::pow(dv2coll_proton_basic_at_mid_y_err, 2.0) +
		   std::pow(dv2MF_proton_basic_at_mid_y_err, 2.0) );
      // dv3total/dt in 4pi
      const double dv3total_in_4pi =
	dv3coll_proton_basic_in_4pi + dv3MF_proton_basic_in_4pi;
      const double dv3total_in_4pi_err =
	std::sqrt( std::pow(dv3coll_proton_basic_in_4pi_err, 2.0) +
		   std::pow(dv3MF_proton_basic_in_4pi_err, 2.0) );
      // dv3total/dt at mid y
      const double dv3total_at_mid_y =
	dv3coll_proton_basic_at_mid_y + dv3MF_proton_basic_at_mid_y;
      const double dv3total_at_mid_y_err =
	std::sqrt( std::pow(dv3coll_proton_basic_at_mid_y_err, 2.0) +
		   std::pow(dv3MF_proton_basic_at_mid_y_err, 2.0) );
      // Integrate to obtain int dvnMF (and enable int dvntotal) in 4pi
      int_dv1MF_proton_basic_in_4pi += dv1MF_proton_basic_in_4pi;
      int_dv2MF_proton_basic_in_4pi += dv2MF_proton_basic_in_4pi;
      int_dv3MF_proton_basic_in_4pi += dv3MF_proton_basic_in_4pi;
      // Integrate to obtain int dvnMF (and enable int dvntotal) at mid y
      int_dv1MF_proton_basic_at_mid_y += dv1MF_proton_basic_at_mid_y;
      int_dv2MF_proton_basic_at_mid_y += dv2MF_proton_basic_at_mid_y;
      int_dv3MF_proton_basic_at_mid_y += dv3MF_proton_basic_at_mid_y;

      /////////////////////////////////////
      // Record dvnMF and int dvnMF, and dvntotal and int dvntotal:
      // dv1MF/dt in 4pi
      g_dv1MF_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv1MF_proton_basic_in_4pi/dt);
      g_dv1MF_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv1MF_proton_basic_in_4pi_err/dt);
      // dv1MF/dt at mid y
      g_dv1MF_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv1MF_proton_basic_at_mid_y/dt);
      g_dv1MF_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv1MF_proton_basic_at_mid_y_err/dt);
      // dv2MF/dt in 4pi
      g_dv2MF_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv2MF_proton_basic_in_4pi/dt);
      g_dv2MF_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv2MF_proton_basic_in_4pi_err/dt);
      // dv2MF/dt at mid y
      g_dv2MF_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv2MF_proton_basic_at_mid_y/dt);
      g_dv2MF_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv2MF_proton_basic_at_mid_y_err/dt);
      // dv3MF/dt in 4pi
      g_dv3MF_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv3MF_proton_basic_in_4pi/dt);
      g_dv3MF_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv3MF_proton_basic_in_4pi_err/dt);
      // dv3MF/dt at mid y
      g_dv3MF_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv3MF_proton_basic_at_mid_y/dt);
      g_dv3MF_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv3MF_proton_basic_at_mid_y_err/dt);
      // dv1total/dt in 4pi
      g_dv1total_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv1total_in_4pi/dt);
      g_dv1total_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv1total_in_4pi_err/dt);
      // dv1total/dt at mid y
      g_dv1total_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv1total_at_mid_y/dt);
      g_dv1total_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv1total_at_mid_y_err/dt);
      // dv2total/dt in 4pi
      g_dv2total_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv2total_in_4pi/dt);
      g_dv2total_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv2total_in_4pi_err/dt);
      // dv2total/dt at mid y
      g_dv2total_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv2total_at_mid_y/dt);
      g_dv2total_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv2total_at_mid_y_err/dt);
      // dv3total/dt in 4pi
      g_dv3total_dt_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, dv3total_in_4pi/dt);
      g_dv3total_dt_proton_basic_integrated_in_4pi_->
	SetPointError(i_MF_entry, 0.0, dv3total_in_4pi_err/dt);
      // dv3total/dt at mid y
      g_dv3total_dt_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, dv3total_at_mid_y/dt);
      g_dv3total_dt_proton_basic_integrated_at_mid_y_->
	SetPointError(i_MF_entry, 0.0, dv3total_at_mid_y_err/dt);
      // int dvnMF in 4pi
      g_v1MF_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, int_dv1MF_proton_basic_in_4pi);	
      g_v2MF_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, int_dv2MF_proton_basic_in_4pi);	
      g_v3MF_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t, int_dv3MF_proton_basic_in_4pi);
      // int dvnMF at mid y
      g_v1MF_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, int_dv1MF_proton_basic_at_mid_y);	
      g_v2MF_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, int_dv2MF_proton_basic_at_mid_y);	
      g_v3MF_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t, int_dv3MF_proton_basic_at_mid_y);
      // int dvntotal in 4pi
      g_v1total_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t,
		 int_dv1coll_proton_basic_in_4pi + int_dv1MF_proton_basic_in_4pi);
      g_v2total_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t,
		 int_dv2coll_proton_basic_in_4pi + int_dv2MF_proton_basic_in_4pi);
      g_v3total_proton_basic_integrated_in_4pi_->
	SetPoint(i_MF_entry, current_t,
		 int_dv3coll_proton_basic_in_4pi + int_dv3MF_proton_basic_in_4pi);
      // int dvntotal at mid y
      g_v1total_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t,
		 int_dv1coll_proton_basic_at_mid_y + int_dv1MF_proton_basic_at_mid_y);
      g_v2total_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t,
		 int_dv2coll_proton_basic_at_mid_y + int_dv2MF_proton_basic_at_mid_y);
      g_v3total_proton_basic_integrated_at_mid_y_->
	SetPoint(i_MF_entry, current_t,
		 int_dv3coll_proton_basic_at_mid_y + int_dv3MF_proton_basic_at_mid_y);
    }
    
  }



  ////////////////////////////////////////////////////////////////////////////////////////
  // Plot and save data on dvn/dt

  const int real_event_equivalent =
    ROOT_file->n_events() * ROOT_file->n_test() * ROOT_file->n_ensembles();

  const int every_nth_entry = 5;

  ///////////////////////////////////////////
  // Plot observables in 4 pi:
  plot_and_save_flow_time_evolution_data
    (g_dv1coll_dt_proton_basic_integrated_in_4pi_,
     g_dv1MF_dt_proton_basic_integrated_in_4pi_,
     g_dv1total_dt_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.15, 0.15,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_dv2coll_dt_proton_basic_integrated_in_4pi_,
     g_dv2MF_dt_proton_basic_integrated_in_4pi_,
     g_dv2total_dt_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.05, 0.05,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_dv3coll_dt_proton_basic_integrated_in_4pi_,
     g_dv3MF_dt_proton_basic_integrated_in_4pi_,
     g_dv3total_dt_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.075, 0.075,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v1coll_proton_basic_integrated_in_4pi_,
     g_v1MF_proton_basic_integrated_in_4pi_,
     g_v1total_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.15, 0.15,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v2coll_proton_basic_integrated_in_4pi_,
     g_v2MF_proton_basic_integrated_in_4pi_,
     g_v2total_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.1, 0.1,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v3coll_proton_basic_integrated_in_4pi_,
     g_v3MF_proton_basic_integrated_in_4pi_,
     g_v3total_proton_basic_integrated_in_4pi_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.075, 0.075,
     every_nth_entry, cfg, true);

  ///////////////////////////////////////////
  // Plot observables at mid y
    plot_and_save_flow_time_evolution_data
    (g_dv1coll_dt_proton_basic_integrated_at_mid_y_,
     g_dv1MF_dt_proton_basic_integrated_at_mid_y_,
     g_dv1total_dt_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.15, 0.15,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_dv2coll_dt_proton_basic_integrated_at_mid_y_,
     g_dv2MF_dt_proton_basic_integrated_at_mid_y_,
     g_dv2total_dt_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.05, 0.05,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_dv3coll_dt_proton_basic_integrated_at_mid_y_,
     g_dv3MF_dt_proton_basic_integrated_at_mid_y_,
     g_dv3total_dt_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.075, 0.075,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v1coll_proton_basic_integrated_at_mid_y_,
     g_v1MF_proton_basic_integrated_at_mid_y_,
     g_v1total_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.15, 0.15,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v2coll_proton_basic_integrated_at_mid_y_,
     g_v2MF_proton_basic_integrated_at_mid_y_,
     g_v2total_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.1, 0.1,
     every_nth_entry, cfg, true);

  plot_and_save_flow_time_evolution_data
    (g_v3coll_proton_basic_integrated_at_mid_y_,
     g_v3MF_proton_basic_integrated_at_mid_y_,
     g_v3total_proton_basic_integrated_at_mid_y_,
     ROOT_file->n_events(), real_event_equivalent,
     ROOT_file->current_t_steps()[0], ROOT_file->max_time(),
     // these values are just guesses
     -0.075, 0.075,
     every_nth_entry, cfg, true);

  

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Finished basic flow time evolution analysis in 4pi and at mid y :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n\n" << color::RESET << std::endl;
}



