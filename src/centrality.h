/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains basic functions used 
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CENTRALITY_H
#define SMASH_ROOT_ANALYSIS_CENTRALITY_H



//////////////////////////////////////////////////////////////////////////////////////////
// A useful class storing values defining a centrality class: lower percentage, upper
// percentage, lower Ncharged multiplicity, upper Ncharged multiplicity
//////////////////////////////////////////////////////////////////////////////////////////

class Centrality {
 public:
  // Constructor I
  Centrality(double percentage_min,
	     double percentage_max,
	     int Ncharged_min,
	     int Ncharged_max)
    // Explicit initialization of const members
    : percentage_min_(percentage_min),
      percentage_max_(percentage_max),
      Ncharged_min_(Ncharged_min),
      Ncharged_max_(Ncharged_max) {
    // Nothing to do here
  }
  // Constructor II
  Centrality(double percentage_min,
	     double percentage_max,
	     int Ncharged_min,
	     int Ncharged_max,
	     double impact_b_mean,
	     double impact_b_5_percent,
	     double impact_b_95_percent)
    // Explicit initialization of const members
    : percentage_min_(percentage_min),
      percentage_max_(percentage_max),
      Ncharged_min_(Ncharged_min),
      Ncharged_max_(Ncharged_max) {
    // Initialize impact parameter values
    impact_b_mean_ = impact_b_mean;
    impact_b_5_percent_ = impact_b_5_percent;
    impact_b_95_percent_ = impact_b_95_percent;
  }
  // Destructor
  virtual ~Centrality() {}

  ///////////////////////////////////////////
  // Member functions
  void set_impact_b_mean(double b_mean) { impact_b_mean_ = b_mean; }
  void set_impact_b_5_percent(double b_5_perc) { impact_b_5_percent_ = b_5_perc; }
  void set_impact_b_95_percent(double b_95_perc) { impact_b_95_percent_ = b_95_perc; }
  

  
  ///////////////////////////////////////////
  // Return functions
  double percentage_min() const { return percentage_min_; }
  double percentage_max() const { return percentage_max_; }
  int Ncharged_min() const { return Ncharged_min_; }
  int Ncharged_max() const { return Ncharged_max_; }
  double impact_b_mean() const { return impact_b_mean_; }
  double impact_b_5_percent() const { return impact_b_5_percent_; }
  double impact_b_95_percent() const { return impact_b_95_percent_; }

  

 private:
  const double percentage_min_;
  const double percentage_max_;
  const int Ncharged_min_;
  const int Ncharged_max_;
  // Impact parameter distribution: the mean and the values at 5% and 95%
  double impact_b_mean_{};
  double impact_b_5_percent_{};
  double impact_b_95_percent_{};
};



#endif
