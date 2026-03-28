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



// Parsing helper
static bool to_bool(const std::string& s) {
  return (s == "true" || s == "1");
}



// Function handling the config entries 
void Config::load(const std::string& config_filename) {

  auto cfg = read_key_value_file(config_filename);

  ///////////////////////////////////////////
  // Directories
  if (cfg.count("Start_directory")) {
    start_directory = std::stoi(cfg["Start_directory"]);
  }

  if (cfg.count("Number_of_directories")) {
    number_of_directories = std::stoi(cfg["Number_of_directories"]);
  }

  ///////////////////////////////////////////
  // Multiplicity
  if (cfg.count("Multiplicity")) {
    multiplicity = to_bool(cfg["Multiplicity"]);
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

  ///////////////////////////////////////////
  // Flow
  if (cfg.count("Flow_basic")) {
    flow_basic = to_bool(cfg["Flow_basic"]);
  }

}



