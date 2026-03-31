/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

#include "read_Particles.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <vector>

#include <TChainElement.h>

#include "constants.h"



//////////////////////////////////////////////////////////////////////////////////////////
// Constructor and destructor implementations
//////////////////////////////////////////////////////////////////////////////////////////

// Default constructor
ReadParticles::ReadParticles() {
  // Initialize fChain to zero
  fChain = 0;

  // The name of the TChain is the same as that of a single SMASH TTree
  TChain *my_chain = new TChain("particles");
  
  // Assume the file is located in the same folder and read the Tree; note: this assumes
  // that the file is called "Particles.root", which is naturally the case for SMASH files
  std::string ROOT_file_address = "./Particles.root";
  my_chain->Add( ROOT_file_address.c_str() );
  // Add the file path to the vector of file paths
  chain_file_paths_.push_back(ROOT_file_address);
  // Establish the metadata file name
  metadata_file_name_ = ROOT_file_address + ".meta";

  // pass the TChain, based on this fChain is initialized
  Init(my_chain);
}



// Alternative constructor 1
ReadParticles::ReadParticles(Config cfg) {
  std::cout << "Forming a TChain from "
	    << cfg.number_of_directories - cfg.start_directory << " files..."
	    << std::endl;
  // Initialize fChain to zero
  fChain = 0;

  // The name of the TChain is the same as that of a single SMASH TTree
  TChain *my_chain = new TChain("particles");

  // Loop over Particles.root files which are assumed to be located in directories named
  // with consecutive integers (as is the case in SMASH); the range of the loop is given
  // by the parameters of the constructor
  for (int i_folder = 0; i_folder < cfg.number_of_directories; i_folder++) {
    const int folder_id = cfg.start_directory + i_folder;

    // provide the address of the file
    std::string ROOT_file_address(Target_Directory);
    ROOT_file_address += "../data/";
    std::stringstream data_1;
    data_1 << (folder_id);
    ROOT_file_address += data_1.str();
    ROOT_file_address += "/Particles.root";

    // print out a message to confirm that files are read correctly
    if (cfg.verbose) {
      std::cout << "Attaching a ROOT file from folder " << ROOT_file_address << std::endl;
    }

    my_chain->Add( ROOT_file_address.c_str() );
    // Add the file path to the vector of file paths
    chain_file_paths_.push_back(ROOT_file_address);
  }

  // Establish the metadata file name
  std::ostringstream meta;
  meta << Target_Directory << "../data/Particles_"
       << cfg.number_of_directories << "_directories" 
       << ".meta";
  metadata_file_name_ = meta.str();

  // pass the TChain, based on this fChain is initialized
  Init(my_chain);  
}



// Alternative constructor 2
ReadParticles::ReadParticles(std::string ROOT_file_address) {
  // Initialize fChain to zero
  fChain = 0;
  
  // The name of the TChain is the same as that of a single SMASH TTree
  TChain *my_chain = new TChain("particles");
  // Use the provided address
  my_chain->Add( ROOT_file_address.c_str() );
  // Add the file path to the vector of file paths
  chain_file_paths_.push_back(ROOT_file_address);

  // Establish the metadata file name
  metadata_file_name_ = ROOT_file_address + ".meta";

  // pass the TChain, based on this fChain is initialized
  Init(my_chain);
}



// Destructor
ReadParticles::~ReadParticles() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}



//////////////////////////////////////////////////////////////////////////////////////////
// Functions to obtain data from the TTree
//////////////////////////////////////////////////////////////////////////////////////////

Int_t ReadParticles::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;

  Int_t return_int = fChain->GetEntry(entry);
  ////////////////////////////////////////
  // Check that the entry does not contain more particles than Max_Particles (which
  // set the size of the single-particle arrays like p0, t, x, pdgcode, etc.)
  if ( npart > Max_Particles ) {
    throw std::runtime_error("\n\n***************************************************\n"
			     "The ROOT entry contains more particles than the size \n"
			     "of used particle arrays.\n"
			     "Increase Max_Particles in constants.h"
			     "\n*************************************************\n\n");
  }
  
  
  return return_int;
}



Long64_t ReadParticles::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}



//////////////////////////////////////////////////////////////////////////////////////////
// The Init() function, called when the selector needs to initialize a new tree or chain.
// Typically, here the branch addresses and branch pointers of the tree will be set. It is
// normally not necessary to make changes to the generated code, but the routine can be
// extended by the user if needed. Init() will be called many times when running on PROOF
// (once per file to be processed).
//////////////////////////////////////////////////////////////////////////////////////////
void ReadParticles::Init(TChain *tree) {
  // Set branch addresses and branch pointers
  if (!tree) return;
  
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ev", &ev, &b_ev);
  fChain->SetBranchAddress("ens", &ens, &b_ens);
  fChain->SetBranchAddress("tcounter", &tcounter, &b_tcounter);
  fChain->SetBranchAddress("npart", &npart, &b_npart);
  fChain->SetBranchAddress("test_p", &test_p, &b_test_p);
  fChain->SetBranchAddress("modus_l", &modus_l, &b_modus_l);
  fChain->SetBranchAddress("current_t", &current_t, &b_current_t);
  fChain->SetBranchAddress("impact_b", &impact_b, &b_impact_b);
  fChain->SetBranchAddress("empty_event", &empty_event, &b_empty_event);
  // It's necessary to put a safeguard around the id branch which is non-standard
  if ( fChain->GetBranch("id") ) {
    fChain->SetBranchAddress("id", &id, &b_id);
    id_branch_exists = true;
  }
  fChain->SetBranchAddress("pdgcode", pdgcode, &b_pdgcode);
  fChain->SetBranchAddress("charge", charge, &b_charge);
  // It's necessary to put safeguards around the formation_time and time_last_collision
  // branches which are non-standard
  if ( fChain->GetBranch("formation_time") ) {
    fChain->SetBranchAddress("formation_time", formation_time, &b_formation_time);
    formation_time_branch_exists = true;
  }
  if ( fChain->GetBranch("time_last_collision") ) {
    fChain->SetBranchAddress
      ("time_last_collision", time_last_collision, &b_time_last_collision);
    time_last_collision_branch_exists = true;
  }
  fChain->SetBranchAddress("p0", p0, &b_p0);
  fChain->SetBranchAddress("px", px, &b_px);
  fChain->SetBranchAddress("py", py, &b_py);
  fChain->SetBranchAddress("pz", pz, &b_pz);
  fChain->SetBranchAddress("t", t, &b_t);
  fChain->SetBranchAddress("x", x, &b_x);
  fChain->SetBranchAddress("y", y, &b_y);
  fChain->SetBranchAddress("z", z, &b_z);
  fChain->SetBranchAddress("E_kinetic_tot", &E_kinetic_tot, &b_E_kinetic_tot);
  fChain->SetBranchAddress("E_fields_tot", &E_fields_tot, &b_E_fields_tot);
  fChain->SetBranchAddress("E_tot", &E_tot, &b_E_tot);
  
  Notify();
}



//////////////////////////////////////////////////////////////////////////////////////////
// The Notify() function is called when a new file is opened. This can be either for a new
// TTree in a TChain or when a new TTree is started when using PROOF. It is normally not
// necessary to make changes to the generated code, but the routine can be extended by the
// user if needed. The return value is currently not used.
//////////////////////////////////////////////////////////////////////////////////////////
Bool_t ReadParticles::Notify() {
  return kTRUE;
}



void ReadParticles::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  
  fChain->Show(entry);
}



void ReadParticles::print_ROOT_file_properties() const {

  std::cout << "\n\n*****************************************************************"
	    << "\n*****************************************************************"
	    << "\nInformation from the ROOT file:"
	    << "\n                   n_entries_ = " << n_entries()
	    << "\n                    n_events_ = " << n_events()
	    << "\n             n_event_outputs_ = " << n_event_outputs()
    	    << "\n          n_event_time_steps_ = " << n_event_time_steps()
	    << "\n"
    	    << "\n             output_interval_ = " << output_interval()
	    << "\n                    max_time_ = " << max_time()
	    << "\n"
	    << "\n                      n_test_ = " << n_test()
	    << "\n                 n_ensembles_ = " << n_ensembles()
	    << "\n"
	    << "\n                  box_length_ = " << box_length()
	    << "\n number_of_particles_at_init_ = " << number_of_particles_at_init()
	    << "\n" << std::endl;

  for (int i = 0; i < static_cast<int>(time_steps().size()); i++) {
    std::cout << "      time_steps[" << i << "] = " << time_steps()[i]
	      << "\t\t current_t_steps[" << i << "] = " << current_t_steps()[i]
	      << std::endl;
  }  

  std::cout << std::endl; 
}



//////////////////////////////////////////////////////////////////////////////////////////
// These functions enable skipping get_properties(), which can be extremely expensive, if
// it already has been computed for the SAME data set.
//////////////////////////////////////////////////////////////////////////////////////////

// Compute the hash (unique fingerprint) of the used TChain
size_t ReadParticles::compute_chain_hash() const {
  TObjArray* files = fChain->GetListOfFiles();
  if (!files) return 0;

  std::string combined;

  for (int i = 0; i < files->GetEntries(); i++) {
    TChainElement* el = (TChainElement*)files->At(i);
    std::string path = el->GetTitle();

    struct stat st_file_info{};
    Long64_t size = 0;
    time_t mtime = 0;
    // Use stat() to fill st_file_info with metadata about the file
    if (stat(path.c_str(), &st_file_info) == 0) {
      size = st_file_info.st_size;
      mtime = st_file_info.st_mtime;
    }

    combined += path;
    combined += std::to_string(size);
    combined += std::to_string(mtime);
  }

  return std::hash<std::string>{}(combined);
}



// Load metadata (if valid) from disk
bool ReadParticles::load_metadata_if_valid(size_t hash, Config cfg) {
  std::ifstream metadata_file(metadata_file_name_);
  if (!metadata_file) return false;

  size_t stored_hash = 0;
  int stored_n_files = 0;

  // We know the structure of the metadata file: the first line is the hash, followed by
  // the number of files used to create the metadata file
  metadata_file >> stored_hash;
  metadata_file >> stored_n_files;

  // Sanity check on chain size
  TObjArray* files = fChain->GetListOfFiles();
  int n_files = files ? files->GetEntries() : 0;

  if (stored_hash != hash || stored_n_files != n_files) {
    return false;
  }

  // Since metadata is valid, load ROOT file properties
  std::string key;
  while (metadata_file >> key) {
    if (key == "n_entries") {
      metadata_file >> ROOT_file_properties_.n_entries_;
    }
    else if (key == "n_events") {
      metadata_file >> ROOT_file_properties_.n_events_;
    }
    else if (key == "n_event_outputs") {
      metadata_file >> ROOT_file_properties_.n_event_outputs_;
    }
    else if (key == "n_event_time_steps") {
      metadata_file >> ROOT_file_properties_.n_event_time_steps_;
    }
    else if (key == "output_interval") {
      metadata_file >> ROOT_file_properties_.output_interval_;
    }
    else if (key == "max_time") {
      metadata_file >> ROOT_file_properties_.max_time_;
    }
    else if (key == "n_test") {
      metadata_file >> ROOT_file_properties_.n_test_;
    }
    else if (key == "n_ensembles") {
      metadata_file >> ROOT_file_properties_.n_ensembles_;
    }
    else if (key == "box_length") {
      metadata_file >> ROOT_file_properties_.box_length_;
    }
    else if (key == "number_of_particles_at_init") {
      metadata_file >> ROOT_file_properties_.number_of_particles_at_init_;
    }
    else if (key == "TIME_STEPS") {
      ROOT_file_properties_.time_steps_.clear();
      // Read the entire rest of the line
      std::string line;
      std::getline(metadata_file, line);
      // Read in elements of the line
      std::istringstream iss(line);
      double x;
      while (iss >> x)
        ROOT_file_properties_.time_steps_.push_back(x);
    }
    else if (key == "CURRENT_T_STEPS") {
      ROOT_file_properties_.current_t_steps_.clear();
      // Read the entire rest of the line
      std::string line;
      std::getline(metadata_file, line);
      // Read in elements of the line
      std::istringstream iss(line);
      double x;
      while (iss >> x)
        ROOT_file_properties_.current_t_steps_.push_back(x);
    }
  }

  if (cfg.verbose) {
    print_ROOT_file_properties();
  }

  return true;
}



// Save the metadata
void ReadParticles::write_metadata(size_t hash) const {
  std::ofstream metadata_file(metadata_file_name_);

  TObjArray* files = fChain->GetListOfFiles();
  int n_files = files ? files->GetEntries() : 0;

  metadata_file << hash << "\n";
  metadata_file << n_files << "\n";

  metadata_file << "n_entries " << ROOT_file_properties_.n_entries_ << "\n";
  metadata_file << "n_events " << ROOT_file_properties_.n_events_ << "\n";
  metadata_file << "n_event_outputs " << ROOT_file_properties_.n_event_outputs_ << "\n";
  metadata_file << "n_event_time_steps "
		<< ROOT_file_properties_.n_event_time_steps_ << "\n";
  metadata_file << "output_interval " << ROOT_file_properties_.output_interval_ << "\n";
  metadata_file << "max_time " << ROOT_file_properties_.max_time_ << "\n";
  metadata_file << "n_test " << ROOT_file_properties_.n_test_ << "\n";
  metadata_file << "n_ensembles " << ROOT_file_properties_.n_ensembles_ << "\n";
  metadata_file << "box_length " << ROOT_file_properties_.box_length_ << "\n";
  metadata_file << "number_of_particles_at_init "
		<< ROOT_file_properties_.number_of_particles_at_init_ << "\n";
  
  metadata_file << "TIME_STEPS ";
  for (double x : ROOT_file_properties_.time_steps_) {
    metadata_file << x << " ";
  }
  metadata_file << "\n";

  metadata_file << "CURRENT_T_STEPS ";
  for (double x : ROOT_file_properties_.current_t_steps_) {
    metadata_file << x << " ";
  }
  metadata_file << "\n";

  std::cout << "\n*********************"
	    << "\nMetadata file created"
	    << "\n*********************" << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////////
// This function updates the private class members: number of entries, number of events,
// number of event outputs, number of steps in an event, output interval, final time of
// the evolution, number of test particles per particle, number of ensembles, box length,
// number of particles at initialization, an array listing time steps in an event, and an
// array listing current_t in an event.
//////////////////////////////////////////////////////////////////////////////////////////
void ReadParticles::get_properties(Config cfg) {
  // Compute the TChain hash
  const size_t hash = compute_chain_hash();
  if ( load_metadata_if_valid(hash, cfg) ) {
    std::cout << "\nLoaded cached ROOT metadata.\n";
    return;
  }

  std::cout << "\n\nMetadata invalid or missing. Running full scan..." << std::endl;

  std::cout << "\n\nRecovering basic data about the simulation..." << std::endl;

  // Access the TTree to get basic information about the events
  
  // Check whether the TTree loaded
  if (fChain == 0) {
    throw std::runtime_error("Failed to load the TTree somehow (fChain).");
  }

  std::cout << "\nGet the number of entries: \n" << std::endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Get the number of entries = all time steps across all events, including ensembles
  ROOT_file_properties_.set_n_entries( fChain->GetEntries() );
  std::cout << "\nn_entries = " << n_entries() << std::endl;

  std::cout << "\nRead ROOT file info: \n" << std::endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Define auxiliary variables used to obtain various quantities; initialize them to
  // absurd values so that it's clear if something is/is not overwritten
  double output_interval_aux = -100.0;
  // max_time can't be initialized to -100 because of how it's used to get output interval
  double max_time_aux = 0.0;
  double t_first_entry = -100.0;
  double current_t_first_entry = -100.0;
  double impact_b_aux = -100.0;
  int ens_first_entry = -100;
  int ensemble_aux = -100;

  // Initialize iterated quantities to zero
  int n_event_steps_aux = 0;
  int n_ensembles_aux = 0;
  std::vector<double> time_steps_aux{};
  std::vector<double> current_t_steps_aux{};

  // Auxiliary bools
  bool this_is_the_first_time_step_of_the_first_event = true;
  bool this_is_the_last_event = false;

  ///////////////////////////////////////////
  std::cout << "\nLoop over all entries \n" << std::endl;

  std::cout << "Loaded entry " << std::endl;
  // Loop over the entries to read of quantities of interest
  for (int i = 0; i < n_entries(); i++) {

    //////////////////////////////////////////////////////////////////////////////////////
    // The structure of this function is dictated by some previous as well as persisting
    // features of the SMASH ROOT output; as a result, sometimes it may look like the code
    // obtains quantities of interest in a roundabout way. It likely indeed does,
    // especially given continuous upgrades to the ROOT output, but we want to maintain
    // backward compatibility.
    //////////////////////////////////////////////////////////////////////////////////////

    // For each entry, load the TTree
    Long64_t i_entry = LoadTree(i);
    if (i_entry < 0) {
      std::cout << "\n\nFailed to load the TTree at i = " << i << std::endl;
      break;
    }
    const Long64_t output_for_fChain_1 = GetEntry(i);

    // Print out entry index
    if (i == 0) {
      std::cout << i;
    } else {
      std::cout << ", " << i ;
    }    

    //////////////////////////////////////////////////////////////////////////////////////
    // Read off some of the quantities at the first time step: box_length_,
    // number_of_particles_at_init_, n_test_, the smallest time, and the impact parameter
    if ( i == 0) {
      // Record the box length
      ROOT_file_properties_.set_box_length(modus_l);
      // We want to know the initial particle number only for the box modus, where we use
      // it to calculate average density
      if ( box_length() > 0 ) {
	ROOT_file_properties_.set_number_of_particles_at_init(npart);
      }
      ROOT_file_properties_.set_n_test(test_p);
      // Record the smallest time; for this, we use the time of the zeroth particle in the
      // first ensemble of the first event (current_t used to have a bug); we also record
      // the current_t.
      t_first_entry = t[0];
      time_steps_aux.push_back(t[0]);
      current_t_first_entry = current_t;
      current_t_steps_aux.push_back(current_t);
      // Record the first impact parameter; this is needed because for modi different than
      // collider, impact_b will always be -1, and the if statement reading off impact_b
      // won't work
      impact_b_aux = impact_b;
      // Record the first ensemble index; this is needed because in the collider mode with
      // Nensembles = 1 and a single used value of the impact parameter, having the same
      // values of t_first_entry and ens_first_entry is the only way to see that the next
      // event has been loaded
      ens_first_entry = ens;
      ensemble_aux = ens;
    }


    
    //////////////////////////////////////////////////////////////////////////////////////
    // Conditions to break the loop

    // Break the loop if another event is loaded: based on the impact parameter.
    // This will work only for the Collider mode (AND if a range of impact parameters is
    // used), so we need to make sure we only invoke this if impact_b_aux is ever updated
    // to a value different than -1.
    if ( impact_b_aux != -1 ) {
      // Is a new value of impact parameter encountered? then break
      if ( impact_b_aux != impact_b ) {
	std::cout << "\n*****************************************************************"
		  << "\nif ( impact_b_aux != impact_b ) breaks the loop"
		  << "\n impact_b_aux = " << impact_b_aux
		  << "\n     impact_b = " << impact_b
		  << "\n\n" << std::endl;
	break;
      }
    }
    
    // Break the loop if another event is loaded: based on current_t.
    // This will work well only for the Box mode, so we need to make sure we only invoke
    // this if box_length is greater than -1.
    // NOTE: this still won't work if the Box modus has Only_Final output.
    if ( box_length() > 0.0 ) {
      // End the loop if an entry from the next event is read in
      if ( (max_time_aux > 0.0) && (current_t <= 0.0) ) {
	std::cout << "\n*****************************************************************"
		  << "\nif ( (max_time_aux > 0.0) && (current_t <= 0.0) ) breaks the loop"
		  << "\n max_time_aux = " << max_time_aux
		  << "\n    current_t = " << current_t
		  << "\n\n" << std::endl;
	break;
      }
    }

    // Break the loop if another event is loaded: based on ens.
    // This breaking condition is needed for Collider mode data with Only_Final option for
    // output and Value option for the impact parameter (instead of a Range).
    // This will work only for SMASH versions 3.2 and higher, as before that release the
    // ens output was flawed.
    // Only apply this if all time steps are the same:
    if ( t[0] == t_first_entry ) {
      // If you hit ens=0, that's the next event now
      if ( (ens == 0) && (ensemble_aux > 0) ) {
	std::cout << "\n*****************************************************************"
		  << "\nif ( (ens == 0) && (ensemble_aux > 0) ) breaks the loop"
		  << "\n            ens = " << ens
		  << "\n   ensemble_aux = " << ensemble_aux
		  << "\n\n" << std::endl;
	break;
      }
    }

    // Break the loop if another event is loaded: based on ens and current_t at the i = 1
    // entry.
    // This breaking condition is needed for Collider mode data with Only_Final option for
    // output and Value option for the impact parameter (instead of a Range) in the case
    // where Nensembles = 1.
    // Check that this is the second entry
    if ( i == 1 ) {
      // Check if there is only one ensemble (otherwise ens would be equal 1)
      if ( ens == 0 ) {
	// Check if there is only one time step
	if ( current_t == current_t_first_entry ) {	
	  std::cout << "\n*******************************************************************"
		    << "\nif ( (i == 1) && (ens == 0) && "
		    << "(current_t == current_t_first_entry) ) breaks the loop"
		    << "\n                     i = " << i
		    << "\n                   ens = " << ens
		    << "\n current_t_first_entry = " << current_t_first_entry
		    << "\n             current_t = " << current_t
		    << "\n\n" << std::endl;
	  break;
	}
      }
    }


    
    // TO DO: Another condition for the Box modus which will work even if Only_Final is used
    
    

    if ( t[0] > t_first_entry ) {
      this_is_the_first_time_step_of_the_first_event = false;
    }

    if ( !this_is_the_first_time_step_of_the_first_event && (t[0] == t_first_entry) ) {
      std::cout << "\n*******************************************************************"
		<< "\nif ( !this_is_the_first_time_step_of_the_first_event && "
		<< "(t[0] == t_first_entry) ) breaks the loop now" << std::endl;
      std::cout << "             t[0] = " << t[0]
		<< "\n\n" << std::endl;
      break;
    }

    
    
    //////////////////////////////////////////////////////////////////////////////////////
    // If the loop continues, read off other quantities

    // Count steps in a single event
    n_event_steps_aux++;

    // Count ensembles in a single event
    // Make sure you only look at the same time step
    if ( t_first_entry == t[0] ) {
      // Make sure you only look at the same event (for OnlyFinal, t_aux = t[0] always)
      if (impact_b_aux == impact_b) {
	n_ensembles_aux++;
      }
    }
    
    // Check that we are beyond entries connected to the first time step of the first
    // event (i.e., n_ensembles does not change anymore)
    if ( n_event_steps_aux > n_ensembles_aux ) {
      // Use the last ensemble to record time
      if ( (n_event_steps_aux % n_ensembles_aux) == 0 ) {
	time_steps_aux.push_back(t[0]);
	current_t_steps_aux.push_back(current_t);
      }
    }    

    // Record impact_b_aux if impact_b is not equal -1 (in the collider mode, this happens
    // at the LAST time step in the output, which is when SMASH updates impact_b from -1
    // to its actual value)
    if ( impact_b != -1 ) {
      impact_b_aux = impact_b;
    }

    // Establish the output interval and the final time of the evolution t_end; we don't
    // use current_t for this, because it had a bug in previous versions of SMASH
    output_interval_aux = std::max( output_interval_aux, t[0] - max_time_aux);
    max_time_aux = std::max( max_time_aux, t[0] );

    // Record ensemble number (works properly only for SMASH 3.2 and higher)
    ensemble_aux = ens;
  } /// for(int i = 0; i < n_entries; i++)


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Check that the output interval makes sense
  if ( output_interval_aux == max_time_aux ) {
    // In this case it's likely that this is an OnlyFinal dataset and output_interval has
    // no meaning
    output_interval_aux = 0;
  }

  
 
  ////////////////////////////////////////////////////////////////////////////////////////
  // Get the number of time steps in a single event; note this takes into the account
  // "Only Final" option in SMASH <---- TO DO: check that this is indeed the case
  /*
  int n_event_steps_aux = 0;
  if ( first_output_interval != output_interval ){
    n_event_steps_aux = std::floor( (max_current_t - 0.000001)/output_interval ) + 2;
  } else {
    n_event_steps_aux = 1;
  }  
  */

  ////////////////////////////////////////////////////////////////////////////////////////
  // Extract information
  ROOT_file_properties_.set_n_ensembles(n_ensembles_aux);
  ROOT_file_properties_.set_n_event_outputs(n_event_steps_aux);
  ROOT_file_properties_.set_n_event_time_steps(n_event_steps_aux / this->n_ensembles());
  ROOT_file_properties_.set_n_events
    (std::floor( this->n_entries()/(this->n_event_time_steps() * this->n_ensembles()) ));

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Update remaining class members
  ////////////////////////////////////////////////////////////////////////////////////////
  ROOT_file_properties_.set_output_interval(output_interval_aux);
  ROOT_file_properties_.set_max_time(max_time_aux);
  ROOT_file_properties_.set_time_steps(time_steps_aux);
  ROOT_file_properties_.set_current_t_steps(current_t_steps_aux);


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Print ROOT file properties
  ////////////////////////////////////////////////////////////////////////////////////////
  if (cfg.verbose) {
    print_ROOT_file_properties();
  }
    
  ////////////////////////////////////////////////////////////////////////////////////////
  // Update the metadata file
  ////////////////////////////////////////////////////////////////////////////////////////
  write_metadata(hash);
}

