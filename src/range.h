/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains basic functions used 
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_RANGE_H
#define SMASH_ROOT_ANALYSIS_RANGE_H



//////////////////////////////////////////////////////////////////////////////////////////
// A class storing doubles defining a range: mix and max, x(5%) and x(95%), etc.
//////////////////////////////////////////////////////////////////////////////////////////

// We template the class so that it can be of any type: int, double, etc.
template <typename T>
class Range {
 public:
  // Default constructor
  Range() = default;
  // Constructor with input
  Range(T min, T max) {
    min_ = min;
    max_ = max;
  }
  // Destructor
  virtual ~Range() {}

  ///////////////////////////////////////////
  // Member functions
  T min() const { return min_; }
  T max() const { return max_; }
  void set_min(T new_min) { min_ = new_min; }
  void set_max(T new_max) { max_ = new_max; }

 private:
  T min_{};
  T max_{};
};



#endif
