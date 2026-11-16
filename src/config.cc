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
#include <stdexcept>
#include <unordered_map>

#include "./range.h"



//////////////////////////////////////////////////////////////////////////////////////////
// Helpers

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



// Parsing helper for vector of ints
static std::vector<int> to_vector_of_ints(const std::string& s) {
  std::vector<int> result;

  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // trim whitespace
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);

    if (!item.empty()) {
      result.push_back(std::stoi(item));
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



// Parsing helper for set of doubles
static std::set<double> to_set_of_doubles(const std::string& s) {
  std::set<double> result;

  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // trim whitespace
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);

    if (!item.empty()) {
      result.insert(std::stod(item));
    }
  }

  return result;
}



// Returns the text between '[' and ']', or throws if the brackets are missing
static std::string inside_brackets(const std::string& s) {
  const size_t open  = s.find('[');
  const size_t close = s.find(']');
  if (open == std::string::npos || close == std::string::npos || close < open) {
    throw std::invalid_argument("Expected a range of the form [min, max], got: " + s);
  }
  return s.substr(open + 1, close - open - 1);
}



// Parsing helper for a range of ints written as [min, max]
static Range<int> to_range_of_ints(const std::string& s) {
  const std::vector<int> values = to_vector_of_ints(inside_brackets(s));
  if (values.size() != 2) {
    throw std::invalid_argument("Expected exactly two values in [min, max], got: " + s);
  }
  if (values[0] > values[1]) {
    throw std::invalid_argument("Range has min > max: " + s);
  }
  return Range<int>(values[0], values[1]);
}



// Parsing helper for a range of doubles written as [min, max]
static Range<double> to_range_of_doubles(const std::string& s) {
  const std::vector<double> values = to_vector_of_doubles(inside_brackets(s));
  if (values.size() != 2) {
    throw std::invalid_argument("Expected exactly two values in [min, max], got: " + s);
  }
  if (values[0] > values[1]) {
    throw std::invalid_argument("Range has min > max: " + s);
  }
  return Range<double>(values[0], values[1]);
}



//////////////////////////////////////////////////////////////////////////////////////////
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
  // Multiplicity
  if (cfg.count("Multiplicity")) {
    multiplicity = to_bool(cfg["Multiplicity"]);
  }
  if (cfg.count("Centrality_class_edges")) {
    centrality_class_edges =
      to_vector_of_doubles(cfg["Centrality_class_edges"]);
  }
  if (cfg.count("Multiplicity_excluded_species")) {
    multiplicity_cuts.pdgs_to_exclude =
      to_unordered_set_of_ints(cfg["Multiplicity_excluded_species"]);
  }
  if (cfg.count("Multiplicity_FXT_frame")) {
    multiplicity_cuts.FXT_frame = to_bool(cfg["Multiplicity_FXT_frame"]);
  }
  // Default cuts
  if (cfg.count("Multiplicity_default_eta_and_pT_cuts")) {
    multiplicity_default_eta_and_pT_cuts =
      to_bool(cfg["Multiplicity_default_eta_and_pT_cuts"]);
    default_multiplicity_eta_and_pT_cuts();
  }
  // Custom cuts
  if (cfg.count("Multiplicity_eta_min")) {
    if (multiplicity_default_eta_and_pT_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_cuts.eta.set_min(std::stod(cfg["Multiplicity_eta_min"]));
  }
  if (cfg.count("Multiplicity_eta_max")) {
    if (multiplicity_default_eta_and_pT_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }
    multiplicity_cuts.eta.set_max(std::stod(cfg["Multiplicity_eta_max"]));
  }
  if (cfg.count("Multiplicity_pT_min")) {
    if (multiplicity_default_eta_and_pT_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_cuts.pT.set_min(std::stod(cfg["Multiplicity_pT_min"]));
  }
  if (cfg.count("Multiplicity_pT_max")) {
    if (multiplicity_default_eta_and_pT_cuts) {
      throw std::runtime_error("Cannot both provide custom multiplicity cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    multiplicity_cuts.pT.set_max(std::stod(cfg["Multiplicity_pT_max"]));
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
  // Time evolution
  if (cfg.count("Flow_only_at_final_output")) {
    flow_only_at_final_output = to_bool(cfg["Flow_only_at_final_output"]);
  }
  if (cfg.count("Flow_output_times")) {
    flow_output_times = to_set_of_doubles(cfg["Flow_output_times"]);
  } 
  // Rapidity coverage
  if (cfg.count("Flow_scale_by_beam_rapidity")) {
    flow_scale_by_beam_rapidity = to_bool(cfg["Flow_scale_by_beam_rapidity"]);
  }
  if (cfg.count("Flow_number_of_rapidity_bins")) {
    flow_number_of_rapidity_bins = std::stoi(cfg["Flow_number_of_rapidity_bins"]);
  }
  if (cfg.count("Flow_y_min")) {
    flow_y_min = std::stod(cfg["Flow_y_min"]);
  }
  if (cfg.count("Flow_y_max")) {
    flow_y_max = std::stod(cfg["Flow_y_max"]);
  }
  if (cfg.count("Flow_y_mid_min")) {
    flow_y_mid_min = std::stod(cfg["Flow_y_mid_min"]);
  }
  if (cfg.count("Flow_y_mid_max")) {
    flow_y_mid_max = std::stod(cfg["Flow_y_mid_max"]);
  }
  // Centrality cuts: a percentage range plus exactly one of an Ncharged or b range
  const bool has_percentage = cfg.count("Flow_centrality_percentage");
  const bool has_Ncharged = cfg.count("Flow_centrality_Ncharged");
  const bool has_impact_b = cfg.count("Flow_centrality_impact_b");

  if (has_Ncharged && has_impact_b) {
    throw std::runtime_error("Give either Flow_centrality_Ncharged or "
			     "Flow_centrality_impact_b, not both.");
  }
  if (has_percentage != (has_Ncharged || has_impact_b)) {
    throw std::runtime_error("Flow_centrality_percentage must be given together with "
			     "Flow_centrality_Ncharged or Flow_centrality_impact_b.");
  }

  if (has_percentage) {
    const Range<double> per = to_range_of_doubles(cfg["Flow_centrality_percentage"]);
    if (per.min() < 0.0 || per.max() > 100.0 || per.min() >= per.max()) {
      throw std::invalid_argument("Invalid Flow_centrality_percentage range.");
    }

    if (has_Ncharged) {
      const Range<int> Nch = to_range_of_ints(cfg["Flow_centrality_Ncharged"]);
      if (Nch.min() < 0) {
	throw std::invalid_argument("Invalid Flow_centrality_Ncharged range.");
      }
      flow_centrality.emplace(per.min(), per.max(), Ncharged_range{Nch.min(), Nch.max()});
    } else {
      const Range<double> b = to_range_of_doubles(cfg["Flow_centrality_impact_b"]);
      if (b.min() < 0.0) {
	throw std::invalid_argument("Invalid Flow_centrality_impact_b range.");
      }
      flow_centrality.emplace(per.min(), per.max(), ImpactB_range{b.min(), b.max()});
    }
  }
  // Default pT cuts
  if (cfg.count("Flow_default_pT_cuts_HADES")) {
    flow_default_pT_cuts_HADES = to_bool(cfg["Flow_default_pT_cuts_HADES"]);
    default_flow_pT_cuts();
  }
  if (cfg.count("Flow_default_pT_cuts_STAR_FXT")) {
    flow_default_pT_cuts_STAR_FXT = to_bool(cfg["Flow_default_pT_cuts_STAR_FXT"]);
    default_flow_pT_cuts();
  }  
  // Custom pT cuts
  if (cfg.count("Flow_proton_pT_min")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    proton_pT_min = std::stod(cfg["Flow_proton_pT_min"]);
  }
  if (cfg.count("Flow_proton_pT_max")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    proton_pT_max = std::stod(cfg["Flow_proton_pT_max"]);
  }
  if (cfg.count("Flow_deuteron_pT_min")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    deuteron_pT_min = std::stod(cfg["Flow_deuteron_pT_min"]);
  }
  if (cfg.count("Flow_deuteron_pT_max")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    deuteron_pT_max = std::stod(cfg["Flow_deuteron_pT_max"]);
  }
  if (cfg.count("Flow_lambda_pT_min")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    lambda_pT_min = std::stod(cfg["Flow_lambda_pT_min"]);
  }
  if (cfg.count("Flow_lambda_pT_max")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    lambda_pT_max = std::stod(cfg["Flow_lambda_pT_max"]);
  }
  if (cfg.count("Flow_pion_pT_min")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    pion_pT_min = std::stod(cfg["Flow_pion_pT_min"]);
  }
  if (cfg.count("Flow_pion_pT_max")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    pion_pT_max = std::stod(cfg["Flow_pion_pT_max"]);
  }
  if (cfg.count("Flow_kaon_pT_min")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    kaon_pT_min = std::stod(cfg["Flow_kaon_pT_min"]);
  }
  if (cfg.count("Flow_kaon_pT_max")) {
    if (flow_default_pT_cuts_HADES || flow_default_pT_cuts_STAR_FXT) {
      throw std::runtime_error("Cannot both provide custom flow pT cuts "
			       "and ask for default cuts.\n"
			       "Adjust the config file.");
    }    
    kaon_pT_max = std::stod(cfg["Flow_kaon_pT_max"]);
  }

}



