/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate yields in the collider mode.
//////////////////////////////////////////////////////////////////////////////////////////

#include "./yields.h"

#include <filesystem>
#include <fstream>
#include <iostream>
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



void Yields::plot_and_save_1D_histogram_wrapper
  (TH1D* h1, const char* x_axis_label, const char* y_axis_label,
   const char* plot_option, const bool put_plots_in_separate_directory)
{
  // If put_plots_in_separate_directory is true, then make sure the directory exists
  if (put_plots_in_separate_directory) {
    // Full path: Target_Directory + directory name
    std::filesystem::path folder_path =
      std::filesystem::path(Target_Directory) / yields_directory_name_;

    std::error_code ec;
    std::filesystem::create_directories(folder_path, ec);
    if (ec) {
      std::cerr << "Warning: could not create directory "
		<< folder_path.string() << ": " << ec.message() << "\n";
    }
  }
  
  // Get the name for the plot and files (log and not log)
  char basic_file_name_h1[Char_Array_Size];
  snprintf(basic_file_name_h1, Char_Array_Size,
	   "%s%s%s_nEvents=%d",
	   Target_Directory,
	   put_plots_in_separate_directory ? yields_directory_name_ : "",
	   h1->GetName(), N_events_);
  char file_name_h1_log[Char_Array_Size];
  snprintf(file_name_h1_log, Char_Array_Size,
	   "%s_log_plot",
	   basic_file_name_h1);

  plot_and_save_1D_histogram(h1,
			     basic_file_name_h1,
			     h1->GetXaxis()->GetXmin(),
			     h1->GetXaxis()->GetXmax(),
			     x_axis_label, y_axis_label, plot_option,
			     false);
}



void Yields::get_dN_dy
  (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo& config_info,
   const Config cfg)
{
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform yields analysis:
  // 1) create dN/dy histograms for chosen particle species
  // 2) identify yields for specific rapidity bin widths around midrapidity
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Starting yields analysis                            :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n"
	    << color::RESET << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Establish dN/dy distribution for chosen particle species
  
  ///////////////////////////////////////////
  // Loop over all entries
  const int progress_message_threshold = 100;
  std::cout << "Print progress message every "
	    << progress_message_threshold << " entries:" << std::endl;

  // Declare variables for a manual check
  int n_protons_0point1{};
  int n_protons_0point2{};
  int n_protons_0point3{};
  int n_protons_0point4{};
  int n_protons_0point5{};

  for (int i_tree_entry = 0; i_tree_entry < ROOT_file->n_entries(); i_tree_entry++) {
    // Print out a message every X entries
    if ( (i_tree_entry % progress_message_threshold) == 0 ) {
      std::cout << "Yields analysis: loading i_tree_entry = "
		<< i_tree_entry << std::endl;
    }

    // Load the TTree data
    Long64_t Long64_t_for_entry = ROOT_file->LoadTree(i_tree_entry);
    if (Long64_t_for_entry < 0) {
      std::cout << "i_tree_entry = " << i_tree_entry << std::endl;
      throw std::runtime_error("Failed to load the TTree at the above entry.");
    }
    long long int l_l_int_for_entry = ROOT_file->GetEntry(i_tree_entry);
    
    // Only calculate yields at the end time of the simulation
    if ( ROOT_file->current_t < config_info.End_Time() ) {
      continue;
    }
    

    
    /////////////////////////////////////////
    // Loop over all particles in the entry
    for (int index = 0; index < ROOT_file->npart; index++) {
      
      ///////////////////////////////////////
      // protons
      if ( ROOT_file->pdgcode[index] == 2212 ) {
	const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > proton_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_proton_dN_dy_.Fill(y_r);

	  // For a manual check
	  if ( std::abs(y_r) < 0.5 ) {
	    n_protons_0point5++;
	    if ( std::abs(y_r) < 0.4 ) {
	      n_protons_0point4++;
	      if ( std::abs(y_r) < 0.3 ) {
		n_protons_0point3++;
		if ( std::abs(y_r) < 0.2 ) {
		  n_protons_0point2++;
		  if ( std::abs(y_r) < 0.1 ) {
		    n_protons_0point1++;
		  }
		}
	      }
	    }
	  }
	}
      }

      ///////////////////////////////////////
      // antiprotons
      if ( ROOT_file->pdgcode[index] == -2212 ) {
	const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > proton_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_antiproton_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // lambdas
      if ( ROOT_file->pdgcode[index] == 3122 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > lambda_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_lambda_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // antilambdas
      if ( ROOT_file->pdgcode[index] == -3122 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > lambda_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_antilambda_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // pi pluses
      if ( ROOT_file->pdgcode[index] == 211 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > pi_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_pi_plus_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // pi minuses
      if ( ROOT_file->pdgcode[index] == -211 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]);  
	// Apply pT cut
	if ( pT > pi_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_pi_minus_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // K pluses
      if ( ROOT_file->pdgcode[index] == 321 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]); 
	// Apply pT cut
	if ( pT > kaon_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_kaon_plus_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // K minuses
      if ( ROOT_file->pdgcode[index] == -321 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]); 
	// Apply pT cut
	if ( pT > kaon_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_kaon_minus_dN_dy_.Fill(y_r);	  
	}
      }

      ///////////////////////////////////////
      // phi mesons
      if ( ROOT_file->pdgcode[index] == 333 ) {
        const double pT = transverse_momentum(ROOT_file->px[index], ROOT_file->py[index]); 
	// Apply pT cut
	if ( pT > phi_pT_min_ ) {
	  // Get rapidity
	  const double y_r = rapidity(ROOT_file->pz[index], ROOT_file->p0[index]);
	  // Fill the histogram
	  h_phi_dN_dy_.Fill(y_r);	  
	}
      }

    } // for (int index = 0; index < npart; index++) {

  } // for (int i_tree_entry = 0; i_tree_entry < ROOT_file->n_entries(); i_tree_entry++) {

  std::cout << "\n\nFinished filling yields histograms \n\n" << std::endl;


  
  ///////////////////////////////////////////
  // Normalize all histograms to show averages per real event equivalent
  // (i.e., averages comparable with the experiment)
  const int real_event_equivalent =
    ROOT_file->n_events() * ROOT_file->n_test() * ROOT_file->n_ensembles();
  h_proton_dN_dy_.Scale(1.0/real_event_equivalent);
  h_antiproton_dN_dy_.Scale(1.0/real_event_equivalent);
  h_lambda_dN_dy_.Scale(1.0/real_event_equivalent);
  h_antilambda_dN_dy_.Scale(1.0/real_event_equivalent);
  h_pi_plus_dN_dy_.Scale(1.0/real_event_equivalent);
  h_pi_minus_dN_dy_.Scale(1.0/real_event_equivalent);
  h_kaon_plus_dN_dy_.Scale(1.0/real_event_equivalent);
  h_kaon_minus_dN_dy_.Scale(1.0/real_event_equivalent);
  h_phi_dN_dy_.Scale(1.0/real_event_equivalent);
  // Do the same with the manual check variables
  const double n_protons_0point1_per_event =
    static_cast<double>(n_protons_0point1) / (1.0 * real_event_equivalent);
  const double n_protons_0point2_per_event =
    static_cast<double>(n_protons_0point2) / (1.0 * real_event_equivalent);
  const double n_protons_0point3_per_event =
    static_cast<double>(n_protons_0point3) / (1.0 * real_event_equivalent);
  const double n_protons_0point4_per_event =
    static_cast<double>(n_protons_0point4) / (1.0 * real_event_equivalent);
  const double n_protons_0point5_per_event =
    static_cast<double>(n_protons_0point5) / (1.0 * real_event_equivalent);

  if (cfg.verbose) {
    // Print out the values from the manual check
    std::cout << "Proton yields:"
	      << "\n|y| < 0.1: " << n_protons_0point1_per_event
	      << "\n|y| < 0.2: " << n_protons_0point2_per_event
	      << "\n|y| < 0.3: " << n_protons_0point3_per_event
	      << "\n|y| < 0.4: " << n_protons_0point4_per_event
	      << "\n|y| < 0.5: " << n_protons_0point5_per_event
	      << std::endl;
  }
  
  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Extract particle yields for various midrapidty bin widths

  // Store the yields in vectors
  std::vector<double> rapidity_bin_widths =
    {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const int n_of_rapidity_bin_widths = rapidity_bin_widths.size();
  // The yields will be averaged, therefore we use doubles
  std::vector<double> N_proton(n_of_rapidity_bin_widths, 0.0); 
  std::vector<double> N_antiproton(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_lambda(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_antilambda(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_pi_plus(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_pi_minus(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_kaon_plus(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_kaon_minus(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_phi(n_of_rapidity_bin_widths, 0.0);
  // We also store the information about the errors
  std::vector<double> N_proton_error(n_of_rapidity_bin_widths, 0.0); 
  std::vector<double> N_antiproton_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_lambda_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_antilambda_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_pi_plus_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_pi_minus_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_kaon_plus_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_kaon_minus_error(n_of_rapidity_bin_widths, 0.0);
  std::vector<double> N_phi_error(n_of_rapidity_bin_widths, 0.0);

  ///////////////////////////////////////////
  // Loop over all histogram bins
  
  // histogram bin indexes start at 1
  for (int i = 1; i < (n_of_h_bins_ + 1); i++) {
    // This is common to all dN/dy histograms by construction
    const double bin_center = h_proton_dN_dy_.GetBinCenter(i);

    // Loop over all rapidity cuts; we loop in reverse order to be able to break when
    // the rapidity bin center first does not fulfill the condition
    for (int j = (n_of_rapidity_bin_widths - 1); j >= 0 ; j--) {
      // Based on the known bin structure, with centers at +-0.05, +-0.15, etc., we can
      // easily apply the cuts
      if ( std::abs(bin_center) < rapidity_bin_widths[j] ) {
	// Add bin content and error in quadrature
	N_proton[j] += h_proton_dN_dy_.GetBinContent(i);
	N_proton_error[j] += std::pow(h_proton_dN_dy_.GetBinError(i), 2.0);
	
	N_antiproton[j] += h_antiproton_dN_dy_.GetBinContent(i);
	N_antiproton_error[j] += std::pow(h_antiproton_dN_dy_.GetBinError(i), 2.0);
	
	N_lambda[j] += h_lambda_dN_dy_.GetBinContent(i);
	N_lambda_error[j] += std::pow(h_lambda_dN_dy_.GetBinError(i), 2.0);
	
	N_antilambda[j] += h_antilambda_dN_dy_.GetBinContent(i);
	N_antilambda_error[j] += std::pow(h_antilambda_dN_dy_.GetBinError(i), 2.0);
	
	N_pi_plus[j] += h_pi_plus_dN_dy_.GetBinContent(i);
	N_pi_plus_error[j] += std::pow(h_pi_plus_dN_dy_.GetBinError(i), 2.0);
	
	N_pi_minus[j] += h_pi_minus_dN_dy_.GetBinContent(i);
	N_pi_minus_error[j] += std::pow(h_pi_minus_dN_dy_.GetBinError(i), 2.0);
	
	N_kaon_plus[j] += h_kaon_plus_dN_dy_.GetBinContent(i);
	N_kaon_plus_error[j] += std::pow(h_kaon_plus_dN_dy_.GetBinError(i), 2.0);
	
	N_kaon_minus[j] += h_kaon_minus_dN_dy_.GetBinContent(i);
	N_kaon_minus_error[j] += std::pow(h_kaon_minus_dN_dy_.GetBinError(i), 2.0);
	
	N_phi[j] += h_phi_dN_dy_.GetBinContent(i);
	N_phi_error[j] += std::pow(h_phi_dN_dy_.GetBinError(i), 2.0);
      } else {
	// this means no further entries will satisfy the condition
	break;
      }
    }
  }

  // Get actual errors by taking a square root of the errors vectors contents
  for (int i = 0; i < n_of_rapidity_bin_widths; i++) {
    N_proton_error[i] = std::sqrt( N_proton_error[i] );
    N_antiproton_error[i] = std::sqrt( N_antiproton_error[i] );
    N_lambda_error[i] = std::sqrt( N_lambda_error[i] );
    N_antilambda_error[i] = std::sqrt( N_antilambda_error[i] );
    N_pi_plus_error[i] = std::sqrt( N_pi_plus_error[i] );
    N_pi_minus_error[i] = std::sqrt( N_pi_minus_error[i] );
    N_kaon_plus_error[i] = std::sqrt( N_kaon_plus_error[i] );
    N_kaon_minus_error[i] = std::sqrt( N_kaon_minus_error[i] );
    N_phi_error[i] = std::sqrt( N_phi_error[i] );
  }



  ///////////////////////////////////////////
  // Sanity checks
  std::vector<double> n_protons_0pointX_per_event = {n_protons_0point1_per_event,
    n_protons_0point2_per_event, n_protons_0point3_per_event,
    n_protons_0point4_per_event, n_protons_0point5_per_event};
  const double epsilon = 1e-9;
  for (int i = 0; i < 5; i++) {
    if ( std::abs(N_proton[i] - n_protons_0pointX_per_event[i]) > epsilon ) {
      std::ostringstream oss;
      oss << "Fatal error: N_proton[" << i << "] = " << N_proton[i] << ", "
	  << "n_protons_0point" << i << "_per_event = " << n_protons_0pointX_per_event[i]
	  << ".\nThese should be equal (tested precision: " << epsilon << ").";
      throw std::runtime_error(oss.str());
    }
  }
  


  
  

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Save data
  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Saving yields data                                  :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << color::RESET << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  if (cfg.verbose) {
    // Print out the values
    std::cout << "\nYields:" << std::endl;
    for (int i = 0; i < n_of_rapidity_bin_widths; i++) {
      printf("|y| < %2.1f:   N(p) = %8.4f+-%8.4f   N(p\u0305) = %8.4f+-%8.4f   "
	     "N(Λ) = %8.4f+-%8.4f   N(Λ\u0305) = %8.4f+-%8.4f   N(π+) = %8.4f+-%8.4f   \n"
	     "\t    N(π-) = %8.4f+-%8.4f   N(K+) = %8.4f+-%8.4f   N(K-) = %8.4f+-%8.4f   "
	     "N(φ) = %8.4f+-%8.4f\n",
	     rapidity_bin_widths[i],
	     N_proton[i], N_proton_error[i], N_antiproton[i], N_antiproton_error[i],
	     N_lambda[i], N_lambda_error[i], N_antilambda[i], N_antilambda_error[i],
	     N_pi_plus[i], N_pi_plus_error[i], N_pi_minus[i], N_pi_minus_error[i],
	     N_kaon_plus[i], N_kaon_plus_error[i], N_kaon_minus[i], N_kaon_minus_error[i],
	     N_phi[i], N_phi_error[i]);     
    }
    std::cout << "\n" << std::endl;
  }


  
  ///////////////////////////////////////////
  // Save yields to a file

  // Add data in different files for each rapidity bin width
  FILE *Yields_data[n_of_rapidity_bin_widths];
  for (int i = 0; i < n_of_rapidity_bin_widths; i++) {
    // Establish the file name
    char yields_file_name[Char_Array_Size];
    snprintf(yields_file_name, Char_Array_Size,
	     "%s00_yields_rapidity_dependence_sqrts=%.2f_Deltay=%2.1f_nEvents=%d.txt",
	     Target_Directory, sqrts_, rapidity_bin_widths[i], N_events_);
    // Open the file
    Yields_data[i] = fopen(yields_file_name, "w");
  
    // Add header
    std::fprintf(Yields_data[i], "# Yields data                \n#\n"
		 "#   proton pT min = %3.2f\n"
		 "#   Lambda pT min = %3.2f\n"
		 "#     pion pT min = %3.2f\n"
		 "#     kaon pT min = %3.2f\n"
		 "#      phi pT min = %3.2f\n#\n",
		 proton_pT_min_, lambda_pT_min_, pi_pT_min_, kaon_pT_min_, phi_pT_min_);
    // Mark which rapidity bin width
    std::fprintf(Yields_data[i], "# |y| < %2.1f:\n#\n", rapidity_bin_widths[i]);
    // Add column labels
    std::fprintf(Yields_data[i],
		 "#fit?   pdg1    pdg2   fdd1n   fddn2    value       error\n#\n");

    // Do not ask to fit particles whose yields are zero
    const int fit_proton = std::abs(N_proton[i]) > 0 ? 1 : 0;
    const int fit_antiproton = std::abs(N_antiproton[i]) > 0 ? 1 : 0;
    const int fit_lambda = std::abs(N_lambda[i]) > 0 ? 1 : 0;
    const int fit_antilambda = std::abs(N_antilambda[i]) > 0 ? 1 : 0;
    const int fit_pi_plus = std::abs(N_pi_plus[i]) > 0 ? 1 : 0;
    const int fit_pi_minus = std::abs(N_pi_minus[i]) > 0 ? 1 : 0;
    const int fit_kaon_plus = std::abs(N_kaon_plus[i]) > 0 ? 1 : 0;
    const int fit_kaon_minus = std::abs(N_kaon_minus[i]) > 0 ? 1 : 0;
    const int fit_phi = std::abs(N_phi[i]) > 0 ? 1 : 0;
      
    // Particle yields for that bin width
    std::fprintf(Yields_data[i],
		 "%d      2212    0      1       0       %8.4f   %8.4f\n"
		 "%d      -2212   0      1       0       %8.4f   %8.4f\n"
		 "%d      3122    0      1       0       %8.4f   %8.4f\n"
		 "%d      -3122   0      1       0       %8.4f   %8.4f\n"
		 "%d      211     0      1       0       %8.4f   %8.4f\n"
		 "%d      -211    0      1       0       %8.4f   %8.4f\n"
		 "%d      321     0      1       0       %8.4f   %8.4f\n"
		 "%d      -321    0      1       0       %8.4f   %8.4f\n"
		 "%d      333     0      1       0       %8.4f   %8.4f\n",
		 fit_proton, N_proton[i], N_proton_error[i],
		 fit_antiproton, N_antiproton[i], N_antiproton_error[i],
		 fit_lambda, N_lambda[i], N_lambda_error[i],
		 fit_antilambda, N_antilambda[i], N_antilambda_error[i],
		 fit_pi_plus, N_pi_plus[i], N_pi_plus_error[i],
		 fit_pi_minus, N_pi_minus[i], N_pi_minus_error[i],
		 fit_kaon_plus, N_kaon_plus[i], N_kaon_plus_error[i],
		 fit_kaon_minus, N_kaon_minus[i],N_kaon_minus_error[i],
		 fit_phi, N_phi[i], N_phi_error[i] );
    std::fprintf(Yields_data[i], "#\n#\n");

    fclose(Yields_data[i]);
  }

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Divide all histograms by the bin width to obtain dN/dy vs. y
  h_proton_dN_dy_.Scale(1.0/y_bin_width_);
  h_antiproton_dN_dy_.Scale(1.0/y_bin_width_);
  h_lambda_dN_dy_.Scale(1.0/y_bin_width_);
  h_antilambda_dN_dy_.Scale(1.0/y_bin_width_);
  h_pi_plus_dN_dy_.Scale(1.0/y_bin_width_);
  h_pi_minus_dN_dy_.Scale(1.0/y_bin_width_);
  h_kaon_plus_dN_dy_.Scale(1.0/y_bin_width_);
  h_kaon_minus_dN_dy_.Scale(1.0/y_bin_width_);
  h_phi_dN_dy_.Scale(1.0/y_bin_width_);


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Plot and save the histograms
  plot_and_save_1D_histogram_wrapper(&h_proton_dN_dy_,
				     "rapidity y", "dN_{p}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_antiproton_dN_dy_,
				     "rapidity y", "dN_{#bar{p}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_lambda_dN_dy_,
				     "rapidity y", "dN_{#Lambda}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_antilambda_dN_dy_,
				     "rapidity y", "dN_{#bar{#Lambda}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_pi_plus_dN_dy_,
				     "rapidity y", "dN_{#pi^{+}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_pi_minus_dN_dy_,
				     "rapidity y", "dN_{#pi^{0}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_kaon_plus_dN_dy_,
				     "rapidity y", "dN_{K^{+}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_kaon_minus_dN_dy_,
				     "rapidity y", "dN_{K^{-}}/dy", "h");
  plot_and_save_1D_histogram_wrapper(&h_phi_dN_dy_,
				     "rapidity y", "dN_{#phi}/dy", "h");


  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Finished yields analysis                            :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n\n"
	    << color::RESET << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  
}



