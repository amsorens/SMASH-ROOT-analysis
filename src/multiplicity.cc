/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class calculates multiplicity and centrality classes.
//////////////////////////////////////////////////////////////////////////////////////////

#include "./multiplicity.h"

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
#include "./fourvector.h"
#include "./plots.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"



// passes std::unique_ptr<TH1D> by reference to avoid copy and update it 
void Multiplicity::project_onto_h_impact_b
  (double Nch_min, double Nch_max, std::unique_ptr<TH1D>& h_b, const char* h_b_name) {
  // Convert the range of multiplicity values to bin numbers
  int bin_x_min = h_Nch_and_b_.GetXaxis()->FindBin( Nch_min );
  int bin_x_max = h_Nch_and_b_.GetXaxis()->FindBin( Nch_max );
  // Set the range on the X axis (multiplicity)
  h_Nch_and_b_.GetXaxis()->SetRange(bin_x_min, bin_x_max);
  // Project onto the Y axis (impact parameter)
  h_b.reset(h_Nch_and_b_.ProjectionY(h_b_name));
  h_b->SetDirectory(nullptr);
  // Reset range afterwards
  h_Nch_and_b_.GetXaxis()->SetRange();
}



void Multiplicity::plot_and_save_1D_histogram_wrapper
  (const std::unique_ptr<TH1D>& h1, const char* x_axis_label, const char* y_axis_label,
   const char* plot_option, const bool put_plots_in_separate_directory)
{  
  // If put_plots_in_separate_directory is true, then make sure the directory exists
  if (put_plots_in_separate_directory) {
    // Full path: Target_Directory + directory name
    std::filesystem::path folder_path =
      std::filesystem::path(Target_Directory) / multiplicty_directory_name_;

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
	   "%s%s%s%s_"
	   "etamin=%3.2f_etamax=%3.2f_pTmin=%3.2f_pTmax=%3.2f_nEventsNonempty=%d",
	   Target_Directory,
	   put_plots_in_separate_directory ? multiplicty_directory_name_ : "",
	   h1->GetName(), FXT_frame_ ? "_FXT_frame" : "",
	   eta_.min(), eta_.max(), pT_.min(), pT_.max(), N_events_nonempty_);
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
  plot_and_save_1D_histogram(h1,
			     file_name_h1_log,
			     h1->GetXaxis()->GetXmin(),
			     h1->GetXaxis()->GetXmax(),
			     x_axis_label, y_axis_label, plot_option,
			     true);
}



void Multiplicity::plot_and_save_2D_histogram_wrapper
  (TH2D h2, const char* x_axis_label, const char* y_axis_label,
   const char* plot_option, const bool put_plots_in_separate_directory)
{
  // If put_plots_in_separate_directory is true, then make sure the directory exists
  if (put_plots_in_separate_directory) {
    // Full path: Target_Directory + directory name
    std::filesystem::path folder_path =
      std::filesystem::path(Target_Directory) / multiplicty_directory_name_;

    std::error_code ec;
    std::filesystem::create_directories(folder_path, ec);
    if (ec) {
      std::cerr << "Warning: could not create directory "
		<< folder_path.string() << ": " << ec.message() << "\n";
    }
  }
  
  // Get the name for the plot and files (log and not log)
  char basic_file_name_h2[Char_Array_Size];
  snprintf(basic_file_name_h2, Char_Array_Size,
	   "%s%s%s%s_"
	   "etamin=%3.2f_etamax=%3.2f_pTmin=%3.2f_pTmax=%3.2f_nEventsNonempty=%d",
	   Target_Directory,
	   put_plots_in_separate_directory ? multiplicty_directory_name_ : "",
	   h2.GetName(), FXT_frame_ ? "_FXT_frame" : "",
	   eta_.min(), eta_.max(), pT_.min(), pT_.max(), N_events_nonempty_);
  char file_name_h2_log[Char_Array_Size];
  snprintf(file_name_h2_log, Char_Array_Size,
	   "%s_log_plot",
	   basic_file_name_h2);

  plot_and_save_2D_histogram(h2,
			     basic_file_name_h2,
			     h2.GetXaxis()->GetXmin(),
			     h2.GetXaxis()->GetXmax(),
			     h2.GetYaxis()->GetXmin(),
			     h2.GetYaxis()->GetXmax(),
			     x_axis_label, y_axis_label, plot_option,
			     false);  
  plot_and_save_2D_histogram(h2,
			     file_name_h2_log,
			     h2.GetXaxis()->GetXmin(),
			     h2.GetXaxis()->GetXmax(),
			     h2.GetYaxis()->GetXmin(),
			     h2.GetYaxis()->GetXmax(),
			     x_axis_label, y_axis_label, plot_option,
			     true);
}



std::string Multiplicity::pdg_exclusion_string() const {
  std::ostringstream oss;
  oss << "{";
  bool first = true;
  for (int pdg : pdgs_to_exclude_) {
    if (!first) {
      oss << ", ";
    }
    first = false;
    oss << pdg;
  }
  oss << "}";
  return oss.str();
}



void Multiplicity::multiplicity_and_centrality
  (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo config_info)
{
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform centrality analysis mimicking the way it is done in the experiment, i.e.:
  // 1) create a multiplicity histogram
  // 2) going from the upper end of the histogram (high multiplicity), start adding
  //    numbers of events in the histogram bins until you arrive at 5% of the total, 10%
  //    of the total, etc.
  // 3) create histograms of impact parameters within each centrality class
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  std::cout << color::BLUE
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Starting multiplicity and centrality analysis       :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n"
	    << color::RESET << std::endl;


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Get the number of entries and events
  const Long64_t n_entries = ROOT_file->n_entries();
  const int n_events = ROOT_file->n_events();
  


  ////////////////////////////////////////////////////////////////////////////////////////
  // Establish multiplicity distribution

  ///////////////////////////////////////////
  // Print out parameters of the analysis
  std::cout << "\nWe exclude the following PDG codes (see SMASH particle list):"
	    << std::endl;
  for (const auto& pdg_code : pdgs_to_exclude_) {
    std::cout << pdg_code << std::endl;
  }
  std::cout << "\nWe want the following centrality classes:"
	    << std::endl;
  for (int i = 0; i < (centrality_class_edges_.size() - 1); i++) {
    if (i == 0) {
      printf("%4.1f-%4.1f%%",
	     100 * centrality_class_edges_[i], 100 * centrality_class_edges_[i+1]);
    } else{
      printf(", %4.1f-%4.1f%%",
	     100 * centrality_class_edges_[i], 100 * centrality_class_edges_[i+1]);
    }
  }
  std::cout << std::endl;

  if (FXT_frame_) {
    std::cout << "The analysis is performed in the FXT frame" << std::endl;
  }
  std::cout << "The following cuts are used:"
	    << "\n eta_min = " << eta_.min()
	    << "\n eta_max = " << eta_.max()
	    << "\n  pT_min = " << pT_.min() << " GeV"
	    << "\n  pT_max = " << pT_.max() << " GeV\n" << std::endl;

  const int progress_message_threshold = 1000;
  std::cout << "Print progress message every "
	    << progress_message_threshold << " entries:" << std::endl;
    
  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Loop over all entries

  // In the ROOT file, there are Nensembles entries for each "real" event; moreover, if
  // more than only the final time step is recorded, there is a separate entry for each
  // time step. Therefore, to compute event-by-event multiplicity etc., we cannot rely on
  // ROOT entries to correctly partition the information. Instead, we use a manual switch
  // which counts how many ensembles we have gathered information from and which triggers
  // recording (real) event-level information when the current ensemble number equals the
  // total number of ensembles as given by ROOT_file->n_ensembles().
  int which_ensemble = 0;
  // variable to count n_charged within one (real) event
  int event_n_charged = 0;
  
  for (int i_tree_entry = 0; i_tree_entry < n_entries; i_tree_entry++) { /// loop A1
    // Print out a message every X entries
    if ( (i_tree_entry % progress_message_threshold) == 0 ) {
      std::cout << "Multiplicity analysis: loading i_tree_entry = "
		<< i_tree_entry << std::endl;
    }

    // Load the TTree data
    Long64_t Long64_t_for_entry = ROOT_file->LoadTree(i_tree_entry);
    if (Long64_t_for_entry < 0) {
      std::cout << "i_tree_entry = " << i_tree_entry << std::endl;
      throw std::runtime_error("Failed to load the TTree at the above entry.");
    }
    long long int l_l_int_for_entry = ROOT_file->GetEntry(i_tree_entry);
    
    // Only calculate multiplicity at the end time of the simulation
    if ( ROOT_file->current_t < config_info.End_Time() ) {
      continue;
    }


    
    /////////////////////////////////////////
    // Loop over all particles in the entry

    // variable to count the number of charged particles at midrapity for this entry
    int n_charged = 0;
    
    for (int index = 0; index < ROOT_file->npart; index++) { /// loop A2
      // Skip excluded particles
      if (pdgs_to_exclude_.count(ROOT_file->pdgcode[index])) {
	continue;
      }

      // Get the particle four-momentum in the collider frame.
      // TO DO: this assumes the collider frame in the output; in principle, SMASH can be
      // run in the FXT frame, which should be taken into account here.
      FourVector p_mu(ROOT_file->p0[index], ROOT_file->px[index],
		      ROOT_file->py[index], ROOT_file->pz[index]);

      if (FXT_frame_) {
	// Boost to the fixed-target frame
	p_mu.boost_fourvector(0.0, 0.0, -v_beam_);
      }
           
      // Calculate pseudorapidity
      double eta = pseudorapidity(p_mu.x1(), p_mu.x2(), p_mu.x3());
      // Apply the pseudorapidity cut
      if ( eta < eta_.max() ) {
	if ( eta > eta_.min() ) {

	  // Apply the pT cut
	  double pT = transverse_momentum(p_mu.x1(), p_mu.x2());
	  if ( pT > pT_.min() ) {
	    if ( pT < pT_.max() ) {

	      // Only count charged particles
	      if ( std::abs(ROOT_file->charge[index]) > 0 ) {
		n_charged++;
	      }
	    }
	  } // pT cut
	}
      } // eta cut
    } // for (int index = 0; index < npart; index++) { /// loop A2
    
    // Sum n_charged across ensembles
    event_n_charged += n_charged;

    // Update the ensemble index
    which_ensemble++;

    // Record info and reset ensemble counter
    if ( which_ensemble == ROOT_file->n_ensembles() ) {
      // Only record this event if it was nonempty
      if (event_n_charged > 0) {
	N_events_nonempty_++;
	N_charged_.emplace_back(event_n_charged);
	b_impact_.emplace_back(ROOT_file->impact_b);

	// Update the max N_charged if applicable
	if (event_n_charged > N_charged_max_) {
	  N_charged_max_ = event_n_charged;
	}
      }      
      // Reset variables
      event_n_charged = 0;
      which_ensemble = 0;
    }
    
  } // for (int i_tree_entry = 0; i_tree_entry < n_entries; i_tree_entry++) { /// loop A1


  
  ///////////////////////////////////////////
  // Sanity check
  if ( b_impact_.size() != N_charged_.size() ) {
    std::cerr << "ERROR: size mismatch: N_charged_.size() = " << N_charged_.size()
	      << ", b_impact_.size() = " << b_impact_.size() << std::endl;
    throw std::runtime_error("Size mismatch between N_charged_ and b_impact_.");
  }
  if ( static_cast<int>(N_charged_.size()) != N_events_nonempty_ ) {
    std::cerr << "ERROR: size mismatch: N_charged_.size() = " << N_charged_.size()
	      << ", N_events_nonempty_ = " << N_events_nonempty_ << std::endl;
    throw std::runtime_error("Size mismatch between N_charged_ and N_events_nonempty_.");
  }

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Fill basic histograms
  
  // Now that we know N_charged_max_, we initialize and fill the h_Nch_and_b_ histogram
  h_Nch_and_b_.SetName("h_Ncharged_and_impact_b");
  h_Nch_and_b_.SetTitle("Ncharged and impact parameter");
  // For N_ch: set N_charged_max_ + 1 bins centered on integers from 0 to N_charged_max_,
  // resulting in one bin for each possible value of multiplicity.
  // For b: set 181 bins centered on 0.0, 0.1, 0.2, ..., 18.0.
  h_Nch_and_b_.SetBins(N_charged_max_ + 1, -0.5, N_charged_max_ + 0.5,
		       181, -0.05, 18.05);

  // Fill the h_Nch_and_b_histogram
  for (int i = 0; i < static_cast<int>( N_charged_.size() ); i++) {
    h_Nch_and_b_.Fill(N_charged_[i], b_impact_[i]);
  }  
  
  // Project the h_Nch_and_b_ histogram onto the N_ch axis (integrating over b);
  // this needs to be done this way as ProjectionX() returns a heap-allocated TH1D*
  h_Nch_.reset(h_Nch_and_b_.ProjectionX("h_Ncharged"));
  h_Nch_->SetDirectory(nullptr);

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Identify centrality classes

  ///////////////////////////////////////////
  // Loop over all h_Nch_ bins to read off multiplicity classes.

  // We want to sum entries (events) in each bin
  int summed_events = 0;

  // Index variable to read out the upper edge of the given wanted centrality class; it
  // starts from 1 as index 0 is the lower edge of the first class
  int c_edge_index = 1;

  // Auxiliary variables for extracting edges of the centrality classes and the
  // corresponding multiplicities; the two first are initialized to values relevant for
  // the most central class (the first to be updated)
  double current_percentage_min = 0.0;
  int current_Ncharged_max = N_charged_max_;
  double current_percentage_max{};
  int current_Ncharged_min{};  

  // We loop in reverse order: from the highest value bin in N_charged_;
  // bin indices start at 1 in ROOT
  const double epsilon = 1e-12; // small double to tackle floating-point rounding
  for (int i = h_Nch_->GetNbinsX(); i > 0; i--) {
    const double bin_center = h_Nch_->GetBinCenter(i);
    const int bin_content = h_Nch_->GetBinContent(i);
    summed_events += bin_content;
      
    // Compute fraction of events summed so far
    const double fraction = (1.0 * summed_events) / (1.0 * N_events_nonempty_);
    
    // If the fraction is larger than one of the wanted centrality class edges, record the
    // actual centrality percentage, the relevant multiplicities, and update c_edge_index
    if ( fraction >= (centrality_class_edges_[c_edge_index] - epsilon) ) {
      // Assign minimum Ncharged and the upper limit of the centrality class
      current_Ncharged_min = bin_center;
      // up to 2 signficant digits
      current_percentage_max = std::round(fraction * 10000) / 100.0;

      // Populate centrality_classes_
      centrality_classes_.emplace_back
	(Centrality(current_percentage_min, current_percentage_max,
		    current_Ncharged_min, current_Ncharged_max));
      
      // Check if we summed all events
      if ( fraction >= (1.0 - epsilon) ) {
	break;
      }
      // If not, assign the new values to auxiliary variables and update the index:
      // the next lower edge of centrality class is the same as the current upper edge
      current_percentage_min = current_percentage_max;
      // the next Ncharged_max is one less the current Ncharged_min
      current_Ncharged_max = current_Ncharged_min - 1;
      c_edge_index++;      
    }

    // Safety check (this should never happen)
    if ( c_edge_index == (n_of_centrality_classes_+1) ) {
      throw std::runtime_error("c_edge_index reaches (n_of_centrality_classes_ + 1)!!!");
    }
  }

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Project onto impact parameter histograms based on centrality classes

  const int n_of_extracted_centrality_classes = centrality_classes_.size();

  // Loop over all actual centrality classes
  for (int i = 0; i < n_of_extracted_centrality_classes; i++) {
    char hname[Char_Array_Size];
    snprintf(hname, Char_Array_Size, "h_b_impact_%4.2f_%4.2f",
	     centrality_classes_[i].percentage_min(),
	     centrality_classes_[i].percentage_max());
    project_onto_h_impact_b(centrality_classes_[i].Ncharged_min(),
			    centrality_classes_[i].Ncharged_max(),
			    h_impact_b_in_cent_classes_[i], hname);
  }

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Go over the impact parameter histograms and get impact parameters marking the onset
  // of the lower 5% (min) and upper 95% (max) of all impact parameters

  for (int i = 0; i < n_of_extracted_centrality_classes; i++) {
    centrality_classes_[i].set_impact_b_mean(h_impact_b_in_cent_classes_[i]->GetMean());
    
    // Sum over all entries, excluding under- and overflow bins
    const int total_n_entries = (h_impact_b_in_cent_classes_[i])->Integral();

    // Loop from the lowest value bin in h_impact_b_in_cent_classes_[i]
    // (bin indices start at 1 in ROOT)
    int summed_b_events = 0;
    for (int j = 1; j < ((h_impact_b_in_cent_classes_[i])->GetNbinsX() + 1); j++) {
      const double bin_center = (h_impact_b_in_cent_classes_[i])->GetBinCenter(j);
      const int bin_content = (h_impact_b_in_cent_classes_[i])->GetBinContent(j);
      summed_b_events += bin_content;     

      if ( summed_b_events < (0.05 + 1e-12) * total_n_entries ) {
	centrality_classes_[i].set_impact_b_5_percent(bin_center);
      }
      if ( summed_b_events < (0.95 + 1e-12) * total_n_entries ) {
	centrality_classes_[i].set_impact_b_95_percent(bin_center);
      } 
    }
  }
  

  
  ///////////////////////////////////////////
  // Define an auxiliary TH1D that is normalized for plots
  std::unique_ptr<TH1D>
    h_Nch_normalized(h_Nch_and_b_.ProjectionX("h_Ncharged_normalized"));
  h_Nch_normalized->SetDirectory(nullptr);
  h_Nch_normalized->Scale(1.0/N_events_nonempty_);



  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Save data
  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Saving multiplicity and centrality data             :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << color::RESET << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Plot and save the histograms

  plot_and_save_2D_histogram_wrapper(h_Nch_and_b_, "N_charged", "b_impact", "h");

  plot_and_save_1D_histogram_wrapper(h_Nch_, "N_charged", "N_events", "h");
  plot_and_save_1D_histogram_wrapper(h_Nch_normalized, "N_charged", "N_events", "h");
  
  // Loop over all the centrality classes
  for (int i = 0; i < n_of_extracted_centrality_classes; i++) {
    plot_and_save_1D_histogram_wrapper(h_impact_b_in_cent_classes_[i],
				       "b_impact", "N_events", "h");
  }

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Put centrality information into a txt file
  
  // Get the name for the plot and files (log and not log)
  char centrality_file_name[Char_Array_Size];
  snprintf(centrality_file_name, Char_Array_Size,
	   "%s%s%s_sqrts=%3.1f_etamin=%3.2f_etamax=%3.2f_pTmin=%3.2f_pTmax=%3.2f_"
	   "nEventsNonempty=%d.txt",
	   Target_Directory,
	   "centrality_classes", FXT_frame_ ? "_FXT_frame" : "",
	   sqrts_, eta_.min(), eta_.max(), pT_.min(), pT_.max(), N_events_nonempty_);

  // Open the file
  FILE* cent_file = fopen(centrality_file_name, "w");
  if (!cent_file) {
    perror("Error opening the centrality file");
    exit(EXIT_FAILURE);
  }

  // Add a header
  fprintf(cent_file,
	  "# Centrality classes based on %d events \n"
	  "# at sqrts = %.2f GeV and impact parameter range = %s fm\n"
	  "# with cuts as follows%s:\n"
	  "# etamin=%3.2f, etamax=%3.2f, pTmin=%3.2f, pTmax=%3.2f \n"
	  "# and excluding particles with PDG codes %s\n"
	  "# leading to %d non-empty events \n"
	  "# \n"
	  "#        class  Nch_min Nch_max     <b>   b(5%%)  b(95%%) \n"
	  "# \n",
	  N_events_, sqrts_, config_info.Range_or_Value(),
	  FXT_frame_ ? " (in the FXT frame)" : "",
	  eta_.min(), eta_.max(), pT_.min(), pT_.max(),
	  pdg_exclusion_string().c_str(), N_events_nonempty_);

  // Add data
  for (int i = 0; i < n_of_extracted_centrality_classes; i++) {
    fprintf(cent_file, "   %5.2f-%5.2f:   %5d   %5d   %5.2f   %5.1f   %5.1f  \n",
	    centrality_classes_[i].percentage_min(),
	    centrality_classes_[i].percentage_max(),
	    centrality_classes_[i].Ncharged_min(), centrality_classes_[i].Ncharged_max(),
	    centrality_classes_[i].impact_b_mean(),
	    centrality_classes_[i].impact_b_5_percent(),
	    centrality_classes_[i].impact_b_95_percent());
  }

  // Close the file 
  fclose(cent_file);
	   

  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n: Finished multiplicity and centrality analysis       :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n\n"
	    << color::RESET << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
   
}



