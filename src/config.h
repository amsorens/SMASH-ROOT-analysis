/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file defines how to use the config file
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_CONFIG_H
#define SMASH_ROOT_ANALYSIS_CONFIG_H

#include <string>



struct Config {
  // Directories
  int start_directory = 0;
  int number_of_directories = 10;

  // Analyses
  bool multiplicity = false;
  bool yields = false;

  bool flow_basic = false;

  // Read from file
  void load(const std::string& filename);
};



#endif
