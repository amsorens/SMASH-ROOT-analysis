/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains basic functions used 
//////////////////////////////////////////////////////////////////////////////////////////

#include "./basic_functions.h"



//////////////////////////////////////////////////////////////////////////////////////////
// Kinematic variables
//////////////////////////////////////////////////////////////////////////////////////////

// fixed-target frame beam energy; mN has a default value specified in constants file
double Ekin_from_sqrts(double sqrts, double mN) {
  return (sqrts * sqrts)/(2.0 * mN) - 2.0 * mN;
}

// center-of-mass frame beam energy; mN has a default value specified in constants file
double sqrts_from_Ekin(double Ekin, double mN) {
  return std::sqrt( 2.0 * mN * ( Ekin + 2.0 * mN ) );
}

// center-of-mass frame beam rapidity; mN has a default value specified in constants file
double y_beam_cm(double sqrts, double mN) {
  const double E = sqrts/2.0;
  const double pz = sqrt( E * E - mN * mN );
  return 0.5 * log( (E + pz)/(E - pz) );
}

// center-of-mass frame beam velocity; mN has a default value specified in constants file
double v_beam_cm(double sqrts, double mN) {
  const double E = sqrts/2.0;
  const double pz = sqrt( E * E - mN * mN );
  return pz / E;
}

// kinetic energy
double Ekinetic(double mass, double px, double py, double pz) {
  return std::sqrt(mass * mass + px * px + py * py + pz * pz);
}

// kinetic energy, overload for momentum magnitude
double Ekinetic(double mass, double p) {
  return std::sqrt(mass * mass + p * p);
}

// mass
double mass(double E, double px, double py, double pz) {
  return std::sqrt(E * E - px * px - py * py - pz * pz);
}

// Mandelstam variable s for two particles
double Mandelstam_s(double E1, double px1, double py1, double pz1,
		    double E2, double px2, double py2, double pz2) {
  return std::pow(E1 + E2, 2.0) - std::pow(px1 + px2, 2.0) - std::pow(py1 + py2, 2.0)
    - std::pow(pz1 + pz2, 2.0);
}

// pT
double transverse_momentum(double p_x, double p_y) {
  return std::sqrt( p_x * p_x + p_y * p_y );
}

// momentum space rapidity y
double rapidity(double p_z, double E) {
  return 0.5 * std::log( (E + p_z)/(E - p_z) );
}

double pseudorapidity(double p_x, double p_y, double p_z) {
  double p = std::sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
  return 0.5 * std::log( (p + p_z)/(p - p_z) );
}

// azimuthal angle
double phi_angle(double p_x, double p_y, double p_T) {
  double phi = std::acos( p_x / p_T );
  // account for the fact that std::acos only returns in the interval [0, pi]
  if ( p_y < 0) {
    phi = 2.0 * M_PI - phi;
  }
  return phi;
}

// x-component of the center-of-mass velocity
double velocity_COM_x(double E1, double px1, double E2, double px2) {
  return (px1 + px2) / (E1 + E2);
}

// y-component of the center-of-mass velocity
double velocity_COM_y(double E1, double py1, double E2, double py2) {
  return (py1 + py2) / (E1 + E2);
}

// z-component of the center-of-mass velocity
double velocity_COM_z(double E1, double pz1, double E2, double pz2) {
  return (pz1 + pz2) / (E1 + E2);
}

// magnitude of particle momentum in the center-of-mass frame (usually denoted p^*, k^*,
// or q^*)
double momentum_COM(double E1, double px1, double py1, double pz1,
		    double E2, double px2, double py2, double pz2) {
  const double s = Mandelstam_s(E1, px1, py1, pz1, E2, px2, py2, pz2);
  const double m1 = mass(E1, px1, py1, pz1);
  const double m2 = mass(E2, px2, py2, pz2);

  const double term1 = s - std::pow(m1 + m2, 2.0);
  const double term2 = s - std::pow(m1 - m2, 2.0);

  return 0.5 * std::sqrt( term1 * term2 / s );
}

// magnitude of relative velocity in the center-of-mass frame
double velocity_rel_COM(double E1, double px1, double py1, double pz1,
			double E2, double px2, double py2, double pz2) {
  const double s = Mandelstam_s(E1, px1, py1, pz1, E2, px2, py2, pz2);
  const double m1 = mass(E1, px1, py1, pz1);
  const double m2 = mass(E2, px2, py2, pz2);
  const double p_COM = momentum_COM(E1, px1, py1, pz1, E2, px2, py2, pz2);
  const double E1_COM = Ekinetic(m1, p_COM);
  const double E2_COM = Ekinetic(m2, p_COM);

  return ( p_COM * std::sqrt(s) ) / ( E1_COM * E2_COM );
}



//////////////////////////////////////////////////////////////////////////////////////////
// Histogram helper functions
//////////////////////////////////////////////////////////////////////////////////////////

std::pair<double, double> bin_value_and_error (TProfile *histogram, int bin_index) {
  const int bin_number = histogram->GetBin(bin_index);
  const double bin_value = histogram->GetBinContent(bin_number);
  const double bin_error = histogram->GetBinError(bin_number);

  return std::make_pair(bin_value, bin_error);
}


// overload of the above for std::unique_ptr<TProfile> (passed by reference to avoid copy)
std::pair<double, double> bin_value_and_error (const std::unique_ptr<TProfile>& histogram,
					       int bin_index) {
  const int bin_number = histogram->GetBin(bin_index);
  const double bin_value = histogram->GetBinContent(bin_number);
  const double bin_error = histogram->GetBinError(bin_number);

  return std::make_pair(bin_value, bin_error);
}



// convenience overload for a TH1D histogram
std::pair<double, double> bin_value_and_error (TH1D *histogram, int bin_index) {
  const int bin_number = histogram->GetBin(bin_index);
  const double bin_value = histogram->GetBinContent(bin_number);
  const double bin_error = histogram->GetBinError(bin_number);

  return std::make_pair(bin_value, bin_error);
}



