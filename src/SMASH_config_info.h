/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

#ifndef SMASH_ROOT_ANALYSIS_SMASH_CONFIG_INFO_H
#define SMASH_ROOT_ANALYSIS_SMASH_CONFIG_INFO_H

#include <stdexcept>
#include <vector>



std::string cout_a_vector(std::vector<int> vector);

std::string cout_a_vector(std::vector<double> vector);

void read_in_a_char_array(std::string line_from_the_file,
			  int line_reading_index,
			  char* char_to_be_read_in);

double read_in_a_double(std::string line_from_the_file,
			int line_reading_index);

// wrapper for ints
int read_in_an_int(std::string line_from_the_file,
		   int line_reading_index);

std::vector<double> read_in_a_vector(std::string line_from_the_file,
				     int line_reading_index);


//////////////////////////////////////////////////////////////////////////////////////////
// This class reads off chosen SMASH config file entries
//////////////////////////////////////////////////////////////////////////////////////////

class SMASHConfigInfo {
 public:
  // Constructor
  SMASHConfigInfo () {
    // constructor instructions
  }


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Function that reads off chosen variables from the SMASH config file
  void read_SMASH_config(int starting_folder_number);


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Assign values to class members
  // General
  void set_Modus(char* Modus) {
    for (int m = 0; m < 64; m++){
      Modus_[m] = Modus[m];
    }
  }
  void set_End_Time(double End_Time) { End_Time_ = End_Time; }
  void set_Nevents(int Nevents) { Nevents_ = Nevents; }
  void set_Ensembles(int Ensembles) { Ensembles_ = Ensembles; }
  void set_Testparticles(int Testparticles) { Testparticles_ = Testparticles; }
  void set_Triangular_Range(double Triangular_Range) {
    Triangular_Range_ = Triangular_Range;
  }
  void set_Gaussian_Sigma(double Gaussian_Sigma) { Gaussian_Sigma_ = Gaussian_Sigma; }
  // Box modus
  void set_Length(double Length) { Length_ = Length; }
  void set_Temperature(double Temperature) { Temperature_ = Temperature; }
  void set_number_of_neutrons(int number_of_neutrons) {
    number_of_neutrons_ = number_of_neutrons;
  }
  void set_number_of_protons(int number_of_protons) {
    number_of_protons_ = number_of_protons;
  }
  // Collider modus
  void set_Sqrtsnn(double Sqrtsnn) { Sqrtsnn_ = Sqrtsnn; }
  void set_E_Kin(double E_Kin) { E_Kin_ = E_Kin; }
  void set_Range_or_Value(char* Range_or_Value) {
    for (int m = 0; m < 64; m++){
      Range_or_Value_[m] = Range_or_Value[m];
    }
  }
  // Potentials
  void set_VDF_Sat_rhoB(double VDF_Sat_rhoB) { VDF_Sat_rhoB_ = VDF_Sat_rhoB; }
  void set_VDF_Powers(std::vector<double> VDF_Powers) {
    // Check that VDF_Powers has the right size
    const int VDF_powers_size = VDF_Powers.size();
    if ( VDF_powers_size > 4 ) {
      throw std::runtime_error("Provided VDF_Powers vector has a wrong size!"
			       "\nCurrently, VDF_Powers can only have up to 4 entries.");
    }
    // If the check passes, continue
    for (int i = 0; i < VDF_powers_size; i++) {
      VDF_Powers_[i] = VDF_Powers[i];
    }
  }
  void set_VDF_Coeffs(std::vector<double> VDF_Coeffs) {
    // Check that VDF_Coeffs has the right size
    const int VDF_coeffs_size = VDF_Coeffs.size();
    if ( VDF_coeffs_size > 4 ) {
      throw std::runtime_error("Provided VDF_Coeffs vector has a wrong size!"
			       "\nCurrently, VDF_Coeffs can only have up to 4 entries.");
    }
    // If the check passes, continue
    for (int i = 0; i < VDF_coeffs_size; i++) {
      VDF_Coeffs_[i] = VDF_Coeffs[i];
    }
  }
  // Lattice
  void set_Cell_Number(std::vector<int> Cell_Number) {
    // Check that Cell_Number has the right size
    if ( Cell_Number.size() != 3 ) {
      throw std::runtime_error("Provided Cell_Number vector has a wrong size!"
			       "\nCurrently, this needs to be a 3D vector.");
    }
    // If the check passes, continue
    for (int i = 0; i < 3; i++) {
      Cell_Number_[i] = Cell_Number[i];
    }
  }

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for class members
  // General
  char* Modus() { return Modus_; }
  double End_Time() { return End_Time_; } // in fm
  int Nevents() { return Nevents_; }
  int Ensembles() { return Ensembles_; }
  int Testparticles() { return Testparticles_; }
  double Triangular_Range() {return Triangular_Range_; } // in fm
  double Gaussian_Sigma() { return Gaussian_Sigma_; } // in fm
  // Box modus
  double Length() { return Length_; } // in fm
  double Temperature() { return Temperature_; } // in GeV
  int number_of_neutrons() { return number_of_neutrons_; }
  int number_of_protons() { return number_of_protons_; }
  // Collider modus
  double Sqrtsnn() { return Sqrtsnn_; } // in GeV
  double E_Kin() { return E_Kin_; } // in GeV
  char* Range_or_Value() { return Range_or_Value_; } // in fm
  // Potentials
  double VDF_Sat_rhoB() { return VDF_Sat_rhoB_; }
  std::vector<double> VDF_Powers() { return VDF_Powers_; }
  std::vector<double> VDF_Coeffs() { return VDF_Coeffs_; }
  // Lattice
  std::vector<int> Cell_Number() { return Cell_Number_; }


  
 private:
  ////////////////////////////////////////////////////////////////////////////////////////
  // These properties hold values read off of a SMASH config file.
  // Note that all are given default values of zero, and which ones are overwritten in the
  // process of constructing / modyfing the class object depends on the specific case (for
  // example, if a member Sqrtsnn is nonzero, then member Ekin remains zero).
  ////////////////////////////////////////////////////////////////////////////////////////
  // General
  char Modus_[64];
  double End_Time_ = 0.0;
  int Nevents_ = 0;
  int Ensembles_ = 0;
  int Testparticles_ = 0;
  double Triangular_Range_ = 0.0;
  double Gaussian_Sigma_ = 0.0;
  // Box modus
  double Length_ = 0.0;
  double Temperature_ = 0.0;
  int number_of_neutrons_ = 0;
  int number_of_protons_ = 0;
  // Collider modus
  double Sqrtsnn_ = 0.0;
  double E_Kin_ = 0.0;
  char Range_or_Value_[64];
  // Potentials
  double VDF_Sat_rhoB_ = 0.0;
  // VDF_Powers_ and VDF_Coeffs_ assume there are maximum four entries 
  std::vector<double> VDF_Powers_ = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> VDF_Coeffs_ = {0.0, 0.0, 0.0, 0.0};
  // Lattice
  std::vector<int> Cell_Number_ = {0, 0, 0};
};

#endif
