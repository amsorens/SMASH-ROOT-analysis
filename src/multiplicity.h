/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class calculates multiplicity and centrality classes.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_MULTIPLICITY_H
#define SMASH_ROOT_ANALYSIS_MULTIPLICITY_H

#include <memory>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

#include "./basic_functions.h"
#include "./centrality.h"
#include "./range.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"

#include <iostream>


class Multiplicity {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Default constructor, "implemented" here in the .h file
  // We provide default values of the number of bins and spread in rapidity y for the
  // TProfile histograms, as well as the default value of bool for displaying stats
  Multiplicity(double sqrts,
	       int N_events,
	       std::unordered_set<int> pdgs_to_exclude,
	       std::vector<double> centrality_class_edges,
	       bool FXT_frame,
	       double eta_min,
	       double eta_max,
	       double pT_min,
	       double pT_max)
    // Explicit initialization of const members (terms are initialized in the order of
    // declaration in the class definition, regardless of the order of appearance here):
    : sqrts_(sqrts),
      Ekin_( Ekin_from_sqrts(sqrts_) ),
      y_beam_( y_beam_cm(sqrts_) ),
      v_beam_( v_beam_cm(sqrts_) ),
      N_events_(N_events),
      pdgs_to_exclude_(pdgs_to_exclude),
      centrality_class_edges_(centrality_class_edges),
      n_of_centrality_classes_( centrality_class_edges.size() - 1 ),      
      FXT_frame_(FXT_frame),
      ///////////////////////////////////////
      // Cuts:
      eta_(Range<double>(eta_min, eta_max)),
      pT_(Range<double>(pT_min, pT_max)) {
    //////////////////////////////////////////////////////////////////////////////////////
    // Resize to the number of centrality classes
    h_impact_b_in_cent_classes_.resize(n_of_centrality_classes_);
  }

  // Default destructor
  virtual ~Multiplicity() {}


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Member functions of the class  
  ////////////////////////////////////////////////////////////////////////////////////////

  // Projection of the 2D histogram of multiplicities and impact parameters onto a 1D
  // histogram of impact parameters for given multiplicity classes
  void project_onto_h_impact_b
    (double Nch_min, double Nch_max, std::unique_ptr<TH1D>& h_b, const char* h_b_name);

  void plot_and_save_1D_histogram_wrapper
    (const std::unique_ptr<TH1D>& h1, const char* x_axis_label, const char* y_axis_label,
     const char* plot_option, const bool put_plots_in_separate_directory = true);

  void plot_and_save_2D_histogram_wrapper
    (TH2D h2, const char* x_axis_label, const char* y_axis_label,
     const char* plot_option, const bool put_plots_in_separate_directory = true);

  std::string pdg_exclusion_string() const;
  
  // Function checking for non-empty events and populating h_Nch_multiplicity
  void multiplicity_and_centrality
    (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo config_info);

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for private class members
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // Only define return functions for members used outside of the class
  

  
 private:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Parameters of the simulation and analysis:
  // center-of-mass frame beam energy
  const double sqrts_;
  // fixed-target frame beam energy
  const double Ekin_;
  // center-of-mass frame beam rapidity
  const double y_beam_;
  // center-of-mass frame beam velocity
  const double v_beam_;
  // number of events
  const int N_events_;
  // list of particle pdg codes to exclude from centrality analysis
  const std::unordered_set<int> pdgs_to_exclude_;
  // edges of wanted centrality classes
  const std::vector<double> centrality_class_edges_;
  // number of centrality classes
  const int n_of_centrality_classes_;
  // whether the multiplicity is computed in the FXT frame
  const bool FXT_frame_;
  ///////////////////////////////////////////
  // Cuts:
  // considered rapidity range
  const Range<double> eta_;
  // considered transverse momentum pT range
  const Range<double> pT_;
  ////////////////////////////////////////////////////////////////////////////////////////
  // Event characteristics:
  // number of non-empty events
  int N_events_nonempty_{};
  // number of charged particles between eta_min and eta_max in each nonempty event
  std::vector<int> N_charged_{};
  // the corresponding impact parameter
  std::vector<double> b_impact_{};
  // maximal number of charged particles between eta_min and eta_max across all events
  int N_charged_max_{};
  ////////////////////////////////////////////////////////////////////////////////////////
  // Observables:
  // charged_particle multiplicity and impact parameter distribution
  TH2D h_Nch_and_b_; // this will be initialized during analysis
  // projecton of the h_Nch_and_b_ histogram onto the Nch axis
  std::unique_ptr<TH1D> h_Nch_; // this will be initialized during analysis  
  // impact parameter distribution for a given centrality class
  std::vector<std::unique_ptr<TH1D>> h_impact_b_in_cent_classes_;
  // extracted centrality classes
  std::vector<Centrality> centrality_classes_{};
  ////////////////////////////////////////////////////////////////////////////////////////
  // Handling data:
  // name of optional directory to store some of the results in
  const char multiplicity_directory_name_[Char_Array_Size] = "multiplicity_plots/";
};



#endif
