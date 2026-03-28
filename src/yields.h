/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate yields in the collider mode.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_YIELDS_H
#define SMASH_ROOT_ANALYSIS_YIELDS_H

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

#include "./basic_functions.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"



class Yields {
public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Default constructor
  Yields(double sqrts,
	 int N_events,
	 double proton_pT_min,
	 double lambda_pT_min,
	 double pi_pT_min,
	 double kaon_pT_min,
	 double phi_pT_min)
    // Explicit initialization of const terms (terms are initialized in the order of
    // declaration in the class definition, regardless of the order of appearance here):
    : sqrts_(sqrts),
      Ekin_( Ekin_from_sqrts(sqrts_) ),
      y_beam_( y_beam_cm(sqrts_) ),
      // Each histogram will be constructed to range from -y_range_max_ to y_range_max_,
      // covering at least (-1.5*y_beam_, 1.5*y_beam_), with y_range_max_ being an even
      // multiple of 0.1 (harcoded as y_bin_width_)
      y_range_max_( static_cast<int>(std::ceil(15 * y_beam_)) / 10.0 ),
      // The number of bins for each dN/dy histogram is established such that the
      // histograms, constructed to range from -y_range_max_ to y_range_max_, will all
      // have bins of width 0.1 (hardcoded as y_bin_width_)
      n_of_h_bins_( 2.0 * y_range_max_ / y_bin_width_ ),
      N_events_(N_events),
      ///////////////////////////////////////
      // Cuts:
      proton_pT_min_(proton_pT_min),
      lambda_pT_min_(lambda_pT_min),
      pi_pT_min_(pi_pT_min),
      kaon_pT_min_(kaon_pT_min),
      phi_pT_min_(phi_pT_min),
      // Histograms for observables:
      h_proton_dN_dy_("h_proton_dN_dy", "p dN/dy",
		      n_of_h_bins_, -y_range_max_, y_range_max_),
      h_antiproton_dN_dy_("h_antiproton_dN_dy", "#bar{p} dN/dy",
			  n_of_h_bins_, -y_range_max_, y_range_max_),
      h_lambda_dN_dy_("h_lambda_dN_dy", "#Lambda dN/dy",
		      n_of_h_bins_, -y_range_max_, y_range_max_),
      h_antilambda_dN_dy_("h_antilambda_dN_dy", "#bar{#Lambda} dN/dy",
			  n_of_h_bins_, -y_range_max_, y_range_max_),
      h_pi_plus_dN_dy_("h_pi_plus_dN_dy", "#pi^{+} dN/dy",
		       n_of_h_bins_, -y_range_max_, y_range_max_),
      h_pi_minus_dN_dy_("h_pi_minus_dN_dy", "#pi^{-} dN/dy",
			n_of_h_bins_, -y_range_max_, y_range_max_),
      h_kaon_plus_dN_dy_("h_kaon_plus_dN_dy", "K^{+} dN/dy",
			 n_of_h_bins_, -y_range_max_, y_range_max_),
      h_kaon_minus_dN_dy_("h_kaon_minus_dN_dy", "K^{-} dN/dy",
			  n_of_h_bins_, -y_range_max_, y_range_max_),
      h_phi_dN_dy_("h_phi_dN_dy", "#phi dN/dy",
		   n_of_h_bins_, -y_range_max_, y_range_max_){
    //////////////////////////////////////////////////////////////////////////////////////
    // Show the histogram ranges
    std::cout << "\n      y_beam_ = " << y_beam_
	      << "\n y_bin width_ = " << y_bin_width_
	      << "\n y_range_max_ = " << y_range_max_
	      << "\n n_of_h_bins_ = " << n_of_h_bins_
	      << std::endl;
  }

  // Default destructor
  virtual ~Yields() {}


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Member functions of the class  
  ////////////////////////////////////////////////////////////////////////////////////////

  void plot_and_save_1D_histogram_wrapper
    (TH1D* h1, const char* x_axis_label, const char* y_axis_label,
     const char* plot_option);
  
  void get_dN_dy
    (const std::unique_ptr<ReadParticles>& ROOT_file, SMASHConfigInfo& config_info);

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for private class members
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // Return only those that are used outside of the class
  

  
private:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Parameters of the simulation and analysis:
  // center-of-mass frame beam energy
  const double sqrts_;
  // fixed-target frame beam energy
  const double Ekin_;
  // center-of-mass frame beam rapidity
  const double y_beam_;
  // rapidity bin width, hardcoded to be 0.1
  const double y_bin_width_ = 0.1;
  // histogram range; histograms will cover (-y_range_max_, y_range_max_)
  const double y_range_max_;
  // number of histogram bins, based on y_beam_
  const int n_of_h_bins_;  
  // number of events
  const int N_events_;
  ///////////////////////////////////////////
  // Cuts:
  // considered transverse momentum pT ranges
  const double proton_pT_min_;
  const double lambda_pT_min_;
  const double pi_pT_min_;
  const double kaon_pT_min_;
  const double phi_pT_min_;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  // Observables:
  //
  // Particle yields binned in rapidity, normalized to divide out bin width effects;
  // normalization will only be performed after yields are extracted.
  TH1D h_proton_dN_dy_{};
  TH1D h_antiproton_dN_dy_{};
  TH1D h_lambda_dN_dy_{};
  TH1D h_antilambda_dN_dy_{};
  TH1D h_pi_plus_dN_dy_{};
  TH1D h_pi_minus_dN_dy_{};
  TH1D h_kaon_plus_dN_dy_{};
  TH1D h_kaon_minus_dN_dy_{};
  TH1D h_phi_dN_dy_{};
};

#endif
