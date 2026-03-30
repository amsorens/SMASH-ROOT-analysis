/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file defines how to use the config file
//////////////////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>



// Helper function to read in the config file entries
static std::unordered_map<std::string, std::string>
  read_key_value_file(const std::string& config_filename) {

  std::ifstream config_file(config_filename);
  if (!config_file.is_open()) {
    std::cerr << "Error: could not open config file: " << config_filename << std::endl;
    return {};
  }

  std::unordered_map<std::string, std::string> config_entries;

  std::string line;
  while (std::getline(config_file, line)) {
    // Remove comments
    auto pos_comment = line.find('#');
    if (pos_comment != std::string::npos) {
      line = line.substr(0, pos_comment);
    }
    // Skip empty lines
    if (line.find_first_not_of(" \t") == std::string::npos) {
      continue;
    }
    // Split on ':'
    auto pos = line.find(':');
    if (pos == std::string::npos) {
      continue;
    }

    std::string key = line.substr(0, pos);
    std::string value = line.substr(pos + 1);

    // Trim whitespace
    key.erase(0, key.find_first_not_of(" \t"));
    key.erase(key.find_last_not_of(" \t") + 1);

    value.erase(0, value.find_first_not_of(" \t"));
    value.erase(value.find_last_not_of(" \t") + 1);

    config_entries[key] = value;
  }

  return config_entries;
}



// Parsing helper for bools
static bool to_bool(const std::string& s) {
  return (s == "true" || s == "1");
}



// Parsing helper for unordered set of ints
static std::unordered_set<int> to_unordered_set_of_ints(const std::string& s) {
  std::unordered_set<int> result;

  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // trim whitespace
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);

    if (!item.empty()) {
      result.insert(std::stoi(item));
    }
  }

  return result;
}



// Parsing helper for vector of doubles
static std::vector<double> to_vector_of_doubles(const std::string& s) {
  std::vector<double> result;

  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // trim whitespace
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);

    if (!item.empty()) {
      result.push_back(std::stod(item));
    }
  }

  return result;
}



// Function handling the config entries 
void Config::load(const std::string& config_filename) {

  auto cfg = read_key_value_file(config_filename);

  ///////////////////////////////////////////
  // Verbosity
  if (cfg.count("Verbose")) {
    verbose = to_bool(cfg["Verbose"]);
  }

  

  ///////////////////////////////////////////
  // Directories
  if (cfg.count("Start_directory")) {
    start_directory = std::stoi(cfg["Start_directory"]);
  }

  if (cfg.count("Number_of_directories")) {
    number_of_directories = std::stoi(cfg["Number_of_directories"]);
  }


  
  ///////////////////////////////////////////
  // Flow
  if (cfg.count("Flow_basic")) {
    flow_basic = to_bool(cfg["Flow_basic"]);
  }


  
  ///////////////////////////////////////////
  // Multiplicity
  if (cfg.count("Multiplicity")) {
    multiplicity = to_bool(cfg["Multiplicity"]);
  }
  if (cfg.count("Multiplicity_excluded_species")) {
    multiplicity_excluded_species =
      to_unordered_set_of_ints(cfg["Multiplicity_excluded_species"]);
  }
  if (cfg.count("Centrality_class_edges")) {
    centrality_class_edges =
      to_vector_of_doubles(cfg["Centrality_class_edges"]);
  }  
  if (cfg.count("Multiplicity_FXT_frame")) {
    multiplicity_FXT_frame = to_bool(cfg["Multiplicity_FXT_frame"]);
  }
  if (cfg.count("Multiplicity_default_cuts")) {
    multiplicity_default_cuts = to_bool(cfg["Multiplicity_default_cuts"]);
    default_multiplicity_cuts();
  }
  // Custom cuts
  if (cfg.count("Multiplicity_eta_min")) {
    if (multiplicity_default_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_eta_min = std::stod(cfg["Multiplicity_eta_min"]);
  }
  if (cfg.count("Multiplicity_eta_max")) {
    if (multiplicity_default_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }
    multiplicity_eta_max = std::stod(cfg["Multiplicity_eta_max"]);
  }
  if (cfg.count("Multiplicity_pT_min")) {
    if (multiplicity_default_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_pT_min = std::stod(cfg["Multiplicity_pT_min"]);
  }
  if (cfg.count("Multiplicity_pT_max")) {
    if (multiplicity_default_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_pT_max = std::stod(cfg["Multiplicity_pT_max"]);
  }
  
  
  
  ///////////////////////////////////////////
  // Yields
  if (cfg.count("Yields")) {
    yields = to_bool(cfg["Yields"]);
  }
  if (cfg.count("Yields_proton_pT_min")) {
    yields_proton_pT_min = std::stod(cfg["Yields_proton_pT_min"]);
  }
  if (cfg.count("Yields_lambda_pT_min")) {
    yields_lambda_pT_min = std::stod(cfg["Yields_lambda_pT_min"]);
  }
  if (cfg.count("Yields_pi_pT_min")) {
    yields_pi_pT_min = std::stod(cfg["Yields_pi_pT_min"]);
  }
  if (cfg.count("Yields_kaon_pT_min")) {
    yields_kaon_pT_min = std::stod(cfg["Yields_kaon_pT_min"]);
  }
  if (cfg.count("Yields_phi_pT_min")) {
    yields_phi_pT_min = std::stod(cfg["Yields_phi_pT_min"]);
  }

}



