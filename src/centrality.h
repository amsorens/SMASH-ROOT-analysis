/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains a class for handling centrality
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CENTRALITY_H
#define SMASH_ROOT_ANALYSIS_CENTRALITY_H

#include <stdexcept>

#include "./shared_helpers.h"


//////////////////////////////////////////////////////////////////////////////////////////
// A class storing values defining a centrality class: lower percentage, upper percentage,
// and dependning on the selection: either lower Ncharged multiplicity and upper
// Ncharged multiplicity, or lower impact_b parameter and upper impact_b parameter. For
// Ncharged, the interval is closed [N_min, N_max], such that the next interval should be
// [N_max + 1, ...]. For the impact_b, the interval is half-open, [b_min, b_max).
//////////////////////////////////////////////////////////////////////////////////////////

// Helper structs for choosing whether the class is based on Ncharged or impact_b
struct ImpactB_range {
  double min; // fm
  double max;
};

struct Ncharged_range {
  int min;
  int max;
};



class Centrality {
 public:
  // This enum class helps to reject calls on incorrect centrality selector
  enum class Selection {
    ImpactB,
    Ncharged
  };

  // Constructor I: impact-parameter based
  Centrality(double percentage_min,
	     double percentage_max,
	     ImpactB_range range)
    // Explicit initialization of const members
    : percentage_min_(percentage_min),
      percentage_max_(percentage_max),
      selection_(Selection::ImpactB),
      b_range_(range),
      Nch_range_(Ncharged_range{0, 0}) {
    // Nothing to do here
  }

  // Constructor IIA: charged-particle-number based
  Centrality(double percentage_min,
	     double percentage_max,
	     Ncharged_range range)
    // Explicit initialization of const members
    : percentage_min_(percentage_min),
      percentage_max_(percentage_max),
      selection_(Selection::Ncharged),
      Nch_range_(range),
      b_range_(ImpactB_range{0.0, 0.0}) {
    // Nothing to do here
  }

  // Constructor IIB: charged-particle-number based
  Centrality(double percentage_min,
	     double percentage_max,
	     Ncharged_range range,
	     double impact_b_mean,
	     double impact_b_5_percent,
	     double impact_b_95_percent)
    // Explicit initialization of const members
    : percentage_min_(percentage_min),
      percentage_max_(percentage_max),
      selection_(Selection::Ncharged),
      Nch_range_(range),
      b_range_(ImpactB_range{0.0, 0.0}) {
    // Initialize impact parameter values
    impact_b_mean_ = impact_b_mean;
    impact_b_5_percent_ = impact_b_5_percent;
    impact_b_95_percent_ = impact_b_95_percent;
  }

  // Destructor
  virtual ~Centrality() {}



  ///////////////////////////////////////////
  // Member functions
  bool passes_centrality_cut(double impact_b, int n_charged) const {
    if (selection_ == Selection::ImpactB) {
      return impact_b >= b_range_.min && impact_b < b_range_.max;
    } else if (selection_ == Selection::Ncharged) {
      return n_charged >= Nch_range_.min && n_charged <= Nch_range_.max;
    } else {
      throw std::runtime_error("Unknown centrality selector");
    }
  }

  // Set values
  void set_impact_b_mean(double b_mean) { impact_b_mean_ = b_mean; }
  void set_impact_b_5_percent(double b_5_perc) { impact_b_5_percent_ = b_5_perc; }
  void set_impact_b_95_percent(double b_95_perc) { impact_b_95_percent_ = b_95_perc; }
  

  
  ///////////////////////////////////////////
  // Return functions
  double percentage_min() const { return percentage_min_; }
  double percentage_max() const { return percentage_max_; }

  Selection selection() const { return selection_; }

  double impact_b_min() const {
    if (selection_ != Selection::ImpactB) {
      throw std::logic_error("Centrality::impact_b_min() called on a class with a "
			     "different Selection type");
    }
    return b_range_.min;
  }
  double impact_b_max() const {
    if (selection_ != Selection::ImpactB) {
      throw std::logic_error("Centrality::impact_b_max() called on a class with a "
			     "different Selection type");
    }
    return b_range_.max;
  }
  int n_charged_min() const {
    if (selection_ != Selection::Ncharged) {
      throw std::logic_error("Centrality::n_charged_min() called on a class with a "
			     "different Selection type");
    }
    return Nch_range_.min;
  }
  int n_charged_max() const {
    if (selection_ != Selection::Ncharged) {
      throw std::logic_error("Centrality::n_charged_max() called on a class with a "
			     "different Selection type");
    }
    return Nch_range_.max;
  }

  double impact_b_mean() const { return impact_b_mean_; }
  double impact_b_5_percent() const { return impact_b_5_percent_; }
  double impact_b_95_percent() const { return impact_b_95_percent_; }

  

 private:
  const double percentage_min_;
  const double percentage_max_;
  const Selection selection_;
  const ImpactB_range b_range_;
  const Ncharged_range Nch_range_;

  // Impact parameter distribution: the mean and the values at 5% and 95%
  double impact_b_mean_{};
  double impact_b_5_percent_{};
  double impact_b_95_percent_{};
};



#endif
