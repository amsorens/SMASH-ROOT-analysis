/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This class reads off chosen SMASH config file entries
//////////////////////////////////////////////////////////////////////////////////////////

#include "./SMASH_config_info.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <vector>

#include "./basic_functions.h"



std::string cout_a_vector(std::vector<int> vector) {
  std::stringstream ss;
  const int vector_size = vector.size();
  ss << "[";
  for (int i = 0; i < vector_size; i++) {
    ss << vector[i];
    if (i < (vector_size - 1)) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

std::string cout_a_vector(std::vector<double> vector) {
  std::stringstream ss;
  const int vector_size = vector.size();
  ss << "[";
  for (int i = 0; i < vector_size; i++) {
    ss << vector[i];
    if (i < (vector_size - 1)) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}


void read_in_a_char_array(std::string line_from_the_file,
			  int line_reading_index,
			  char* char_to_be_read_in) {
  // a variable to read in a character from the line
  char character_from_the_line;
  
  // declare an empty char array
  char value_from_the_line[64];
  for (int g = 0; g < 64; g++) {
    value_from_the_line[g] = '\0';
  }

  int value_reading_index = line_reading_index + 1;
  int value_recording_index = 0;
  bool continue_reading_in_the_value = true;

  while ( continue_reading_in_the_value ) {
    character_from_the_line = line_from_the_file.at(value_reading_index);

    if ( (character_from_the_line != ' ') ) {
      value_from_the_line[value_recording_index] = character_from_the_line;
      value_recording_index++;
    }
    
    value_reading_index++;
	      
    if ( value_reading_index == static_cast<int>(line_from_the_file.size()) ) {

      for (int k = 0; k < 64; k++){		
        char_to_be_read_in[k] = value_from_the_line[k];
      }
      continue_reading_in_the_value = false;		
    }
	      
  } /// while (continue_reading_in_the_value)
}


double read_in_a_double(std::string line_from_the_file,
			int line_reading_index) {
  double double_to_be_read_in = 0.0;

  // a variable to read in a character from the line
  char character_from_the_line;
  
  // declare an empty char array
  char value_from_the_line[64];
  for (int g = 0; g < 64; g++){
    value_from_the_line[g] = '\0';
  }

  int value_reading_index = line_reading_index + 1;
  int value_recording_index = 0;
  bool continue_reading_in_the_value = true;

  while ( continue_reading_in_the_value ) {
    character_from_the_line = line_from_the_file.at(value_reading_index);

    if ( (character_from_the_line != ' ') ) {
      value_from_the_line[value_recording_index] = character_from_the_line;
      value_recording_index++;
    }
	    
    value_reading_index++;
	      
    if ( value_reading_index == static_cast<int>(line_from_the_file.size()) ) {
      double_to_be_read_in = atof(value_from_the_line);
      continue_reading_in_the_value = false;		
    }
	      
  } /// while (continue_reading_in_the_value)

  return double_to_be_read_in;
}


// wrapper for ints
int read_in_an_int(std::string line_from_the_file,
		   int line_reading_index) {
  double aux_double_to_be_read_in = read_in_a_double(line_from_the_file, line_reading_index);
  return static_cast<int>(aux_double_to_be_read_in);
}


std::vector<double> read_in_a_vector(std::string line_from_the_file,
				     int line_reading_index) {
  std::vector<double> vector_to_be_read_in;
  vector_to_be_read_in.resize(0);

  // a variable to read in a character from the line
  char character_from_the_line;
  
  // declare an empty char array
  char value_from_the_line[64];
  for (int i = 0; i < 64; i++){
    value_from_the_line[i] = '\0';
  }

  int value_reading_index = line_reading_index + 1;
  int value_recording_index = 0;
  bool continue_reading_in_the_value = true;

  int count_recorded_entries = 0;

  while ( continue_reading_in_the_value ) {
    character_from_the_line = line_from_the_file.at(value_reading_index);	    
    value_reading_index++;

    if ( character_from_the_line == '[' ) {
      continue;
    }

    if ( character_from_the_line == ' ' ) {
      continue;
    }

    // If neither of these is true, one needs to continue reading the value
    if ( (character_from_the_line != ',') && (character_from_the_line != ']') ) {
      value_from_the_line[value_recording_index] = character_from_the_line;
      value_recording_index++;
    }

    // If either of these is true, one needs to record the read in value
    if ( (character_from_the_line == ',') || (character_from_the_line == ']') ) {
      vector_to_be_read_in.push_back( atof(value_from_the_line) );
      count_recorded_entries++;
      
      // Reset auxiliary variables
      value_recording_index = 0;
      for (int i = 0; i < 64; i++){
	value_from_the_line[i] = '\0';
      }		
    }

    // Stop reading if you reach the end of the vector
    if ( character_from_the_line == ']' ){
      continue_reading_in_the_value = false;
    }
	      
  } /// while (continue_reading_in_the_value)

  return vector_to_be_read_in;
}



void SMASHConfigInfo::read_SMASH_config(int starting_folder_number) {

  ////////////////////////////////////////////////////////////////////////////////////////
  // Define auxiliary variables that will hold chosen values read from the config file
  ////////////////////////////////////////////////////////////////////////////////////////

  // In order of appearance in the config file; we initialize all the numbers to -100.
  char Modus[64];
  double End_Time = -100.0;
  int Nevents = -100;
  int Ensembles = -100;
  int Testparticles = -100;
  double Triangular_Range = -100.0;
  double Gaussian_Sigma = - 100.0;
  double Length = -100.0;
  double Temperature = -100.0;
  int number_of_neutrons = -100;
  int number_of_protons = -100;
  double Sqrtsnn = -100.0;
  double E_Kin = -100.0;
  char Range_or_Value[64];
  double VDF_Sat_rhoB = -100.0;
  //std::vector<double> VDF_Powers = {-100.0, -100.0, -100.0, -100.0};
  //std::vector<double> VDF_Coeffs = {-100.0, -100.0, -100.0, -100.0};
  std::vector<int> Cell_Number = {-100, -100, -100};
  // set the char arrays to empty
  for (int m = 0; m < 64; m++) {
    Modus[m] = '\0';
    Range_or_Value[m] = '\0';
  }



  ////////////////////////////////////////////////////////////////////////////////////////
  // Open the config file
  ////////////////////////////////////////////////////////////////////////////////////////

  // Declare a string object with an appropriate name; here we assume a certain structure
  // of folders containig analysis, data, etc.
  std::string config_file_address(Target_Directory);
  config_file_address += "../data/";
  std::stringstream folder_no;
  folder_no << starting_folder_number;
  config_file_address += folder_no.str();
  config_file_address += "/config.yaml";
  //std::cout << "\n\nThe file to read in is \n" << config_file_address << std::endl;

  // Read in the file
  std::ifstream config_file( config_file_address.c_str() );

  
  // check whether the config file exists
  bool the_config_file_exists;
  if( !config_file ) { 
    the_config_file_exists = false; // The file was not found.
  } else {                          // If the file was found, then file is non-0.
    the_config_file_exists = true;  // The file was found.
  }
  if ( !the_config_file_exists ) {
    std::stringstream error_message;
    error_message << "\nThe file with the following address doesn't exist:\n"
		  << config_file_address << std::endl;
    throw std::invalid_argument( error_message.str() );
  }


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Define the auxiliary variables for reading the file
  
  // a variable to read in a line from the file
  std::string line_from_the_file;
  // a variable to read in a character from a line
  char character_from_the_line;
  // a variable to hold a formed word
  char word_from_the_line[64];

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Loop over the file entries

  while ( !config_file.eof() ) {
    // get a line from the file
    getline(config_file, line_from_the_file);


    bool continue_reading_the_line = true;

    // discard an empty line (sometimes spurious empty lines are read in)
    if ( line_from_the_file.empty() ) {    
      continue_reading_the_line = false; // this is a throw-away line
    }
    
    // resetting the word_from_the_line to empty
    for (int m = 0; m < 64; m++){
      word_from_the_line[m] = '\0';
    }

    int line_reading_index = 0;
    int word_forming_index = 0;
    bool word_forming_bool = true;
    bool first_valid_character_found = false;

    

    while ( continue_reading_the_line ) {

      character_from_the_line = line_from_the_file.at(line_reading_index);
      //std::cout << character_from_the_line;

      // Avoid some traps
      if ( character_from_the_line == '#' ){
	continue_reading_the_line = false;
	word_forming_bool = false;
	continue;
      }
      if ( character_from_the_line == '{' ){
	continue_reading_the_line = false;
	word_forming_bool = false;
	continue;
      }	
      if ( character_from_the_line == ':' ){
	word_forming_bool = false;
      }


      // Form the word
      if ( word_forming_bool && (character_from_the_line != ' ') ){
	first_valid_character_found = true;
	word_from_the_line[word_forming_index] = character_from_the_line;
	word_forming_index++;	  
      }


	
      if ( !word_forming_bool ) {

        ////////////////////////////////////////////////////////////////////////////////
	// Read in variables of interest

	//std::cout << word_from_the_line << std::endl;
	//std::cin.get();
	
	// Modus
	if ( strncmp(word_from_the_line, "Modus", 5) == 0 ) {
	  read_in_a_char_array(line_from_the_file, line_reading_index, Modus);
	  // assign the value to the class member
	  this->set_Modus(Modus);
	  continue_reading_the_line = false;	  
	}

	// End_Time
	if ( strncmp(word_from_the_line, "End_Time", 8) == 0 ) {	    
	  End_Time = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_End_Time(End_Time);
	  continue_reading_the_line = false;	  
	}

	// Nevents
	if ( strncmp(word_from_the_line, "Nevents", 7) == 0 ) {
	  Nevents = read_in_an_int(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Nevents(Nevents);
	  continue_reading_the_line = false;	
	}

	// Ensembles
	if ( strncmp(word_from_the_line, "Ensembles", 9) == 0 ) {
	  Ensembles = read_in_an_int(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Ensembles(Ensembles);
	  continue_reading_the_line = false;	
	}

	// Testparticles
	if ( strncmp(word_from_the_line, "Testparticles", 13) == 0 ) {
	  Testparticles = read_in_an_int(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Testparticles(Testparticles);
	  continue_reading_the_line = false;	  
	}

	// Triangular_Range
	if ( strncmp(word_from_the_line, "Triangular_Range", 16) == 0 ) {
	  Triangular_Range = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Triangular_Range(Triangular_Range);
	  continue_reading_the_line = false;	  
	}

	// Gaussian_Sigma
	if ( strncmp(word_from_the_line, "Gaussian_Sigma", 14) == 0 ) {
	  Gaussian_Sigma = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Gaussian_Sigma(Gaussian_Sigma);
	  continue_reading_the_line = false;	  
	}

	// Length
	if ( strncmp(word_from_the_line, "Length", 6) == 0 ) {
	  Length = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Length(Length);
	  continue_reading_the_line = false;	
	}

	// Temperature
	if ( strncmp(word_from_the_line, "Temperature", 11) == 0 ) {
	  Temperature = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Temperature(Temperature);
	  continue_reading_the_line = false;	
	}

	// number_of_neutrons
	if ( strncmp(word_from_the_line, "2112", 4) == 0 ) {
	  number_of_neutrons = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_number_of_neutrons(number_of_neutrons);
	  continue_reading_the_line = false;	
	}

	// number_of_protons
	if ( strncmp(word_from_the_line, "2212", 4) == 0 ) {
	  number_of_protons = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_number_of_protons(number_of_protons);
	  continue_reading_the_line = false;	
	}

	// Sqrtsnn
	if ( strncmp(word_from_the_line, "Sqrtsnn", 7) == 0 ) {
	  Sqrtsnn = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_Sqrtsnn(Sqrtsnn);
	  continue_reading_the_line = false;	
	}

	// E_Kin
	if ( strncmp(word_from_the_line, "E_Kin", 5) == 0 ) {
	  E_Kin = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_E_Kin(E_Kin);
	  continue_reading_the_line = false;	
	}

	// (impact parameter) Range
	if ( (strncmp(word_from_the_line, "Range", 5) == 0) ||
	     (strncmp(word_from_the_line, "Value", 5) == 0) ) {
	  read_in_a_char_array(line_from_the_file, line_reading_index, Range_or_Value);
	  // assign the value to the class member
	  this->set_Range_or_Value(Range_or_Value);
	  continue_reading_the_line = false;	
	}

	// (VDF) Sat_rhoB
	if ( strncmp(word_from_the_line, "Sat_rhoB", 8) == 0 ) {
	  VDF_Sat_rhoB = read_in_a_double(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_VDF_Sat_rhoB(VDF_Sat_rhoB);
	  continue_reading_the_line = false;	
	}

	// (VDF) Powers
	if ( strncmp(word_from_the_line, "Powers", 6) == 0 ) {
	  std::vector<double> aux = read_in_a_vector(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_VDF_Powers(aux);
	  continue_reading_the_line = false;	
	}

	// (VDF) Coeffs
	if ( strncmp(word_from_the_line, "Coeffs", 6) == 0 ) {
	  std::vector<double> aux = read_in_a_vector(line_from_the_file, line_reading_index);
	  // assign the value to the class member
	  this->set_VDF_Coeffs(aux);
	  continue_reading_the_line = false;	
	}

	// (lattice) Cell_Number
	if ( strncmp(word_from_the_line, "Cell_Number", 11) == 0 ) {
	  std::vector<double> aux = read_in_a_vector(line_from_the_file, line_reading_index);
	  for (int i = 0; i < 3; i++) {
	    Cell_Number[i] = static_cast<int>(aux[i]); 
	  }
	  // assign the value to the class member
	  this->set_Cell_Number(Cell_Number);
	  continue_reading_the_line = false;	
	}

	

	////////////////////////////////////////////////////////////////////////////////
	// After the last considered option we should make sure that
	continue_reading_the_line = false;	  

      } /// if ( !word_forming_bool )

      line_reading_index++;
      if ( line_reading_index == static_cast<int>(line_from_the_file.size()) ){
	continue_reading_the_line = false;
      }
	  
    } /// while ( continue_reading_the_line )
    
    //std::cout << "\n\n\nWe're moving on to the next line!\n\n\n" << std::endl;
    //std::cin.get();
  } /// while ( !config_file.eof() ) {
  
  config_file.close();
  
  //std::cout << "\n\n\nWe've finished reading the file!\n\n\n" << std::endl;

  // Update kinematic variables; note that these functions use a default value of the
  // nucleon mass mN specified in constants file
  if ( strncmp(Modus_, "Collider", 8) == 0 ) {
    if ( E_Kin_ == 0.0 ) {
      const double E_Kin = Ekin_from_sqrts(Sqrtsnn_);
      this->set_E_Kin(E_Kin);
    }
    if ( Sqrtsnn_ == 0.0 ) {
      const double Sqrtsnn = sqrts_from_Ekin(E_Kin_);
      this->set_Sqrtsnn(Sqrtsnn);
    }
  }

  std::cout << "\n\n*****************************************************************"
	    << "\n*****************************************************************"
	    << "\nInformation from the config file:"
	    << "\n                       Modus_ = " << Modus_
	    << "\n                    End_Time_ = " << End_Time_
	    << "\n                     Nevents_ = " << Nevents_
	    << "\n                   Ensembles_ = " << Ensembles_
	    << "\n               Testparticles_ = " << Testparticles_
	    << "\n            Triangular_Range_ = " << Triangular_Range_
	    << "\n              Gaussian_Sigma_ = " << Gaussian_Sigma_
	    << "\n"
	    << "\n                      Length_ = " << Length_
	    << "\n                 Temperature_ = " << Temperature_
	    << "\n          number_of_neutrons_ = " << number_of_neutrons_
	    << "\n           number_of_protons_ = " << number_of_protons_
	    << "\n"
	    << "\n                     Sqrtsnn_ = " << Sqrtsnn_
	    << "\n                       E_Kin_ = " << E_Kin_
	    << "\n              Range_or_Value_ = " << Range_or_Value_
	    << "\n"
	    << "\n                VDF_Sat_rhoB_ = " << VDF_Sat_rhoB_
	    << "\n                  VDF_Powers_ = " << cout_a_vector(VDF_Powers_)
	    << "\n                  VDF_Coeffs_ = " << cout_a_vector(VDF_Coeffs_)
	    << "\n"
	    << "\n                 Cell_Number_ = " << cout_a_vector(Cell_Number_)
	    << "\n" << std::endl;
  
}

