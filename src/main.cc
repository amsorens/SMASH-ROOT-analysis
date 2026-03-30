/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

#include <chrono>   // for measuring the code execution time with std::chrono
#include <fstream>  // for std::ofstream etc.
#include <iomanip>  // for std::put_time
#include <iostream> // for std::cout etc.
#include <sstream>  // for std::ostringstream

#include "./config.h"
#include "./constants.h"
#include "./multiplicity.h"
#include "./read_Particles.h"
#include "./SMASH_config_info.h"
#include "./yields.h"



int main () {
  std::cout << color::BLUE
	    << "\n\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .."
	    << "\n:                                                     :"
	    << "\n: Welcome to Agnieszka's SMASH ROOT analysis code!    :"
	    << "\n:                                                     :"
	    << "\n: Copyright (c) 2026 Agnieszka Sorensen               :"
	    << "\n:                                                     :"
	    << "\n .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..\n\n"
	    << color::RESET
	    << std::endl;

  ///////////////////////////////////////////
  // Establish what level of verbose statements from ROOT we want
  gErrorIgnoreLevel = kWarning;
  note_msg(std::string("The ROOT warning level used is ") +
	   error_level_to_string(gErrorIgnoreLevel));
  std::cout << std::endl;

  
    
  ///////////////////////////////////////////
  // Start measuring time
  auto t_start = std::chrono::high_resolution_clock::now(); 


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Setup the analysis calculation
  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Read in the config file
  Config cfg;
  cfg.load("./analysis_config.txt");

  

  ///////////////////////////////////////////
  // Initialize a ReadParticles object with the Start_directory and Number_of_directories
  // to attach as specified in the config file.
  // Note: A (unique) pointer is used as the size of the ReadParticles object is expected
  // to be very large.
  std::unique_ptr<ReadParticles> ROOT_file =
    std::make_unique<ReadParticles>(cfg.start_directory, cfg.number_of_directories);

  // Measure the time it takes to initialize
  auto t_init = std::chrono::high_resolution_clock::now();

  

  ///////////////////////////////////////////
  // Get the properties of the ROOT file and SMASH config info
  ROOT_file->get_properties();

  SMASHConfigInfo SMASH_cfg_info;
  SMASH_cfg_info.read_SMASH_config(cfg.start_directory);

  // Measure the time it takes to get properties and info
  auto t_info = std::chrono::high_resolution_clock::now();
  

  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Perform chosen analysis procedures
  ////////////////////////////////////////////////////////////////////////////////////////

  // When the code finishes, check whether any analyses or tests have been performed
  bool any_analysis_performed = false;
  bool any_tests_performed = false;


  std::cout << "\n\n*****************************************************************"
	    << "\n*****************************************************************"
	    << "\nPerforming analyses...\n" << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Basic flow
  if ( cfg.flow_basic ) {
    std::cout << "Basic flow" << std::endl;
    any_analysis_performed = true;
  }


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Multiplicity
  if ( cfg.multiplicity ) {
    if ( strcmp(SMASH_cfg_info.Modus(), "Collider") != 0) {
      throw std::runtime_error("Data is not from a Collider modus \n"
			       "Cannot run multiplicity analysis!");
    }

    Multiplicity multiplicity_analysis(SMASH_cfg_info.Sqrtsnn(),
				       ROOT_file->n_events(),
				       // Get cuts from the config
				       cfg.multiplicity_excluded_species,
				       cfg.centrality_class_edges,
				       cfg.multiplicity_FXT_frame,
				       cfg.multiplicity_eta_min,
				       cfg.multiplicity_eta_max,
				       cfg.multiplicity_pT_min,
				       cfg.multiplicity_pT_max);
    
    multiplicity_analysis.multiplicity_and_centrality(ROOT_file, SMASH_cfg_info);
    any_analysis_performed = true;    
  }


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Yields
  if ( cfg.yields ) {
    if ( strcmp(SMASH_cfg_info.Modus(), "Collider") != 0) {
      throw std::runtime_error("Data is not from a Collider modus \n"
			       "Cannot run yields analysis!");
    }

    Yields yields_analysis(SMASH_cfg_info.Sqrtsnn(),
			   ROOT_file->n_events(),
			   // Get cuts from the config
			   cfg.yields_proton_pT_min, cfg.yields_lambda_pT_min,
			   cfg.yields_pi_pT_min, cfg.yields_kaon_pT_min,
			   cfg.yields_phi_pT_min);
    
    yields_analysis.get_dN_dy(ROOT_file, SMASH_cfg_info, cfg);
    any_analysis_performed = true;
  }



  ////////////////////////////////////////////////////////////////////////////////////////
  // Finish the run
  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Check whether any analyses or tests have been perfomed
  if ( !any_analysis_performed ) {
    if ( !any_tests_performed ) {
      std::cout << "\n\n\nNo analyses or tests were selected to be perfomed."
		<< "\nReview the config file.\n\n\n"
		<< std::endl;
    }
  }



  ////////////////////////////////////////////////////////////////////////////////////////
  // Stop measuring time & calculate duration of each stage
  auto t_end = std::chrono::high_resolution_clock::now();

  auto initialize = std::chrono::duration<double>(t_init - t_start);
  auto get_info = std::chrono::duration<double>(t_info - t_init);
  auto execution = std::chrono::duration<double>(t_end - t_info);
  auto total = std::chrono::duration<double>(t_end - t_start);

  // Output the duration times
  std::cout << "Initialization time: " << initialize.count() << " [s]" << std::endl;
  std::cout << " Retrieve info time: " << get_info.count() << " [s]" << std::endl;
  std::cout << "     Execution time: " << execution.count() << " [s]" << std::endl;
  std::cout << "         Total time: " << total.count() << " [s]" << std::endl;

  // Get current date and time
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm bt = *std::localtime(&now);
  std::ostringstream oss;
  oss << std::put_time(&bt, "%Y-%m-%d %X");

  // Output the duration along with date and time into a file (in append mode)
  std::ofstream output_file("../../000_execution_log.txt", std::ios_base::app);
  if (output_file.is_open()) {
    output_file << "Date and Time: " << oss.str() << std::endl;
    //output_file << "Nevents: " << ROOT_file->n_events() << std::endl;
    output_file << "Initialization time: " << initialize.count() << " [s]" << std::endl;
    output_file << " Retrieve info time: " << get_info.count() << " [s]" << std::endl;
    output_file << "     Execution time: " << execution.count() << " [s] for " << std::endl;
    if (cfg.flow_basic) {
      output_file << "                     basic flow" << std::endl; 
    }
    if (cfg.multiplicity) {
      output_file << "                     multiplicity" << std::endl; 
    }
    if (cfg.yields) {
      output_file << "                     yields" << std::endl; 
    }
    output_file << "         Total time: " << total.count() << " [s]" << std::endl;
    output_file << "---------------------------------------" << std::endl;
    output_file.close();
    std::cout << "Execution time appended to '000_execution_log.txt'" << std::endl;
  } else {
    std::cerr << "Unable to open the execution log!" << std::endl;
  }

  

  return 0;
  
} // main()



