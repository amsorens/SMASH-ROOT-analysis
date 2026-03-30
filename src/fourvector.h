/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class stores components of a fourvector, the corresponding Lorentz invariant, and
// the spatial components of the velocity of the frame in which the fourvector components
// are defined. It also enables performing boosts.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_FOURVECTOR_H
#define SMASH_ROOT_ANALYSIS_FOURVECTOR_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>



class FourVector {
public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Default constructor
  FourVector(double x0, double x1, double x2, double x3,
	     // default frame velocity is for the lab frame
	     double frame_velocity_1 = 0.0,
	     double frame_velocity_2 = 0.0,
	     double frame_velocity_3 = 0.0)
    : x0_(x0),
      x1_(x1),
      x2_(x2),
      x3_(x3),
      Lorentz_invariant_((x0_ * x0_) - (x1_ * x1_) - (x2_ * x2_) - (x3_ * x3_)),
      frame_v1_(frame_velocity_1),
      frame_v2_(frame_velocity_2),
      frame_v3_(frame_velocity_3) {
    // Nothing to do here
  }

  // Copy constructor
  FourVector(const FourVector& other)
    : x0_(other.x0_),
      x1_(other.x1_),
      x2_(other.x2_),
      x3_(other.x3_),
      Lorentz_invariant_(other.Lorentz_invariant_),
      frame_v1_(other.frame_v1_),
      frame_v2_(other.frame_v2_),
      frame_v3_(other.frame_v3_) {
    // Nothing to do here
  }

  // Default destructor
  virtual ~FourVector() {}


  
  double frame_velocity_magnitude() {
    return std::sqrt(frame_v1_ * frame_v1_ + frame_v2_ * frame_v2_ +
		     frame_v3_ * frame_v3_);
  }


  
  void boost_fourvector(double boost_v1, double boost_v2, double boost_v3) {
    // Absolute value of the boost velocity
    const double boost_v =
      std::sqrt(boost_v1 * boost_v1 + boost_v2 * boost_v2 + boost_v3 * boost_v3);
    if ( boost_v >= 1.0 ) {
      throw std::runtime_error("\n\n***************************************************\n"
			       "The provided boost velocity is larger than c"
			       "\n*************************************************\n\n");
    }
    
    // The gamma factor, with a guard against small unstable values
    const double gamma = (boost_v > 1e-9) ? 1.0/std::sqrt(1.0 - boost_v * boost_v) : 1.0;
    // The product \bm{boost_v} \cdot \bm{x}
    const double scalar_product = boost_v1 * x1_ + boost_v2 * x2_ + boost_v3 * x3_;
    // \frac{\gamma - 1}{v^2}
    const double gamma_minus_one_over_v_squared = (gamma - 1.0) / (boost_v * boost_v);
    // The outer bracket in the spatial component transformation
    // \bm{r}' = \bm{r} + (\frac{\gamma - 1}{v^2} (\bm{v} \cdot \bm{r}) - \gamma t) \bm{v}
    const double bracket_term =
      gamma_minus_one_over_v_squared * scalar_product - gamma * x0_;

    // Transform the fourvector; note that we can immediately update the fourvector
    // components as we have already calculated all terms that could have mixed the old
    // and new (boosted) values.
    // NOTE: The order of boosting cannot be changed (x0 first, then xi) or it will break
    // the formulas!!! If needed, introduce temporary variables.
    // Time component transformation t' = \gamma ( t - \bm{v} \cdot \bm{r} )
    x0_ = gamma * (x0_ - scalar_product);
    // Spatial components transformation as \bm{r}' above
    x1_ = x1_ + bracket_term * boost_v1;
    x2_ = x2_ + bracket_term * boost_v2;
    x3_ = x3_ + bracket_term * boost_v3;

    // Calculate the Lorentz invariant and make sure it's preserved
    const double new_Lorentz_inv = (x0_ * x0_) - (x1_ * x1_) - (x2_ * x2_) - (x3_ * x3_);
    if ( std::abs(new_Lorentz_inv - Lorentz_invariant_) > 1e-9 ) {
      std::cout << "\nold_Lorentz_inv = " << Lorentz_invariant_
		<< "\nnew_Lorentz_inv = " << new_Lorentz_inv << "\n" << std::endl;
      throw std::runtime_error("\n\n***************************************************\n"
			       "The new Lorentz invariant is different than the old one!!"
			       "\n*************************************************\n\n");
    }

    // Update frame velocity using 3D relativistic velocity addition:
    // \bm{v}_{old} \cdot \bm{v}_{boost}
    const double v_dot_product =
      frame_v1_ * boost_v1 + frame_v2_ * boost_v2 + frame_v3_ * boost_v3;
    // The outer round bracket in the numerator of
    // \bm{u}' = \frac{ \frac{\bm{u}}{\gamma} +
    //                  (1 + \frac{\gamma}{\gamma + 1} (\bm{v} \cdot \bm{u})) \bm{v}}
    //                {1 + \bm{v} \cdot \bm{u}}
    const double vel_bracket_term = 1.0 + (gamma/(gamma + 1.0)) * v_dot_product;
    // The denominator term in the formula above
    const double denominator = 1.0 / (1.0 + v_dot_product);
    frame_v1_ = (frame_v1_/gamma + vel_bracket_term * boost_v1) * denominator;
    frame_v2_ = (frame_v2_/gamma + vel_bracket_term * boost_v2) * denominator;
    frame_v3_ = (frame_v3_/gamma + vel_bracket_term * boost_v3) * denominator;

    // Check against superluminal frame velocities
    if ( frame_velocity_magnitude() >= 1.0 ) {
      throw std::runtime_error("\n\n***************************************************\n"
			       "The resulting frame velocity is larger than c"
			       "\n*************************************************\n\n");
    }
  }


  
  // Here, the boost velocity is given with respect to the lab frame, but the frame that
  // we boost from (i.e., in which we have the fourvector components) is already moving
  // with some velocity with respect to the lab frame, too. Instead of computing what is
  // the boost velocity in the moving frame, we just perform two boosts.
  void boost_fourvector_from_any_frame(double boost_vx, double boost_vy, double boost_vz,
				       double frame_vx, double frame_vy, double frame_vz) {
    // Apply two boosts: first into the lab frame, then into the new frame
    boost_fourvector(-frame_vx, -frame_vy, -frame_vz);
    boost_fourvector(boost_vx, boost_vy, boost_vz);
  }
  


  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for private class members
  double x0() { return x0_; }
  double x1() { return x1_; }
  double x2() { return x2_; }
  double x3() { return x3_; }
  double frame_velocity_1() { return frame_v1_; }
  double frame_velocity_2() { return frame_v2_; }
  double frame_velocity_3() { return frame_v3_; }
  double Lorentz_invariant() { return Lorentz_invariant_; }
  

  
private:
  // Fourvector components
  double x0_;
  double x1_;
  double x2_;
  double x3_;
  // Lorentz invariant
  const double Lorentz_invariant_;
  // Spatial components of the velocity of the frame in which x0_, x1_, x2_, x3_ are given
  double frame_v1_;
  double frame_v2_;
  double frame_v3_;
};



#endif
