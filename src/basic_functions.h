/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains basic functions used 
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_BASIC_FUNCTIONS_H
#define SMASH_ROOT_ANALYSIS_BASIC_FUNCTIONS_H

#include <TH1.h>
#include <TProfile.h>

#include "./constants.h"



//////////////////////////////////////////////////////////////////////////////////////////
// Messages:
//////////////////////////////////////////////////////////////////////////////////////////
void note_msg(const std::string& msg);

void warning_msg(const std::string& msg);

void success_msg();



//////////////////////////////////////////////////////////////////////////////////////////
// Kinematic variables
//////////////////////////////////////////////////////////////////////////////////////////

// fixed-target frame beam energy
double Ekin_from_sqrts(double sqrts, double mN = Nucleon_Mass);

// center-of-mass frame beam energy
double sqrts_from_Ekin(double Ekin, double mN = Nucleon_Mass);

// center-of-mass frame beam rapidity
double y_beam_cm(double sqrts, double mN = Nucleon_Mass);

// center-of-mass frame beam velocity
double v_beam_cm(double sqrts, double mN = Nucleon_Mass);

// kinetic energy
double Ekinetic(double mass, double px, double py, double pz);

// kinetic energy, overload for momentum magnitude
double Ekinetic(double mass, double p);

// mass
double mass(double E, double px, double py, double pz);

// Mandelstam variable s for two particles
double Mandelstam_s(double E1, double px1, double py1, double pz1,
		    double E2, double px2, double py3, double pz2);

// pT
double transverse_momentum(double p_x, double p_y);

// momentum space rapidity y
double rapidity(double p_z, double E);

double pseudorapidity(double p_x, double p_y, double p_z);

// azimuthal angle
double phi_angle(double p_x, double p_y, double p_T);

// x-component of the center-of-mass velocity
double velocity_COM_x(double E1, double px1, double E2, double px2);

// y-component of the center-of-mass velocity
double velocity_COM_y(double E1, double py1, double E2, double py2);

// z-component of the center-of-mass velocity
double velocity_COM_z(double E1, double pz1, double E2, double pz2);


// magnitude of particle momentum in the center-of-mass frame (usually denoted p^*, k^*,
// or q^*)
double momentum_COM(double E1, double px1, double py1, double pz1,
		    double E2, double px2, double py3, double pz2);


// magnitude of relative velocity in the center-of-mass frame
double velocity_rel_COM(double E1, double px1, double py1, double pz1,
			double E2, double px2, double py3, double pz2);



//////////////////////////////////////////////////////////////////////////////////////////
// Histogram helper functions
//////////////////////////////////////////////////////////////////////////////////////////

std::pair<double, double> bin_value_and_error (TProfile *histogram, int bin_index);

// overload of the above for std::unique_ptr<TProfile> (passed by reference to avoid copy)
std::pair<double, double> bin_value_and_error (const std::unique_ptr<TProfile>& histogram,
					       int bin_index);

// convenience overload for a TH1D histogram
std::pair<double, double> bin_value_and_error (TH1D *histogram, int bin_index);



#endif
