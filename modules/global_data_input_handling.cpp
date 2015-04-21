#include "global_data.h"
#include "global_data_utils.h"

namespace gdu = ::global_data_utils;

namespace {

using gdu::make_quark;
using gdu::quark_check;
using gdu::make_operator_list;

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// simplifies and cleans read_parameters function
void lattice_input_data_handling (const std::string path_output,
    const std::string name_lattice, const std::string path_config, int Lt,
    int Lx, int Ly, int Lz) {
  try{
    if(Lt < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lt\""
          << " is mendatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\n\ttemporal lattice extend .................. " << Lt
        << "\n";
    //
    if(Lx < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lx\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in x direction .... " << Lx
        << "\n";
    //
    if(Ly < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Ly\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in y direction .... " << Ly
        << "\n";
    //
    if(Lz < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lz\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in z direction .... " << Lz
        << "\n\n";
    std::cout << "\tEnsemble ...................................... " <<
      name_lattice << std::endl;
    std::cout << "\tResults will be saved to path:\n\t\t"
        << path_output << "/" << std::endl;
    std::cout << "\tConfigurations will be read from:\n\t\t"
        << path_config << "/" << std::endl;
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

// *****************************************************************************
// simplifies and cleans read_parameters function
void eigenvec_perambulator_input_data_handling (
    const int number_of_eigen_vec, const std::string path_eigenvectors,
    const std::string name_eigenvectors, const std::string path_perambulators,
    const std::string name_perambulators) {

  try{
    if(number_of_eigen_vec < 1){
      std::cout << "\ninput file error:\n" << "\toption \"number_of_eigen_vec\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else{
      std::cout << "\tnumber of eigen vectors .................. "
          << number_of_eigen_vec << "\n";
    }
    std::cout << "\tEigenvectors will be read from files:\n\t\t"
        << path_eigenvectors << "/" << name_eigenvectors
        << "\".eigenvector.t.config\"\n";
    std::cout << "\tPerambulators will be read from files:\n\t\t"
        << path_perambulators << "/" << name_perambulators
        << "\".rnd_vec.scheme.t_sink.config\"\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void quark_input_data_handling (
    const std::vector<std::string> quark_configs, std::vector<quark>& quarks) {
  try{
    // Transform each configured quark into a quark via make_quark, 
    // inserting each object into the quark vector.
    std::transform(quark_configs.begin(), quark_configs.end(),
        std::back_inserter(quarks), make_quark);
    // setting id's in quarks
    size_t quark_counter = 0;
    for(auto&q: quarks){
      q.id = quark_counter;
      quark_counter++;
    }
    // checking the contents for correctness
    std::for_each(quarks.begin(), quarks.end(), quark_check);
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void operator_input_data_handling (
    const std::vector<std::string> operator_list_configs,
    std::vector<Operator_list>& operator_list) {
  try{
    // Transform each configured quark into a quark via make_quark,
    // inserting each object into the quark vector.
    //std::transform(operator_list_configs.begin(), operator_list_configs.end(),
    //    std::back_inserter(operator_list), make_operator_list);
    for(auto op_list : operator_list_configs)
      operator_list.push_back(make_operator_list(op_list));
    // TODO write a check for correctness of input    
  }
  catch(std::exception& e){
    std::cout << "operator_input_data_handling: " << e.what() << "\n";
    exit(0);
  }
}

// *****************************************************************************
/// @brief Makes an operator list object from a string
void correlator_input_data_handling (
    const std::vector<std::string>& correlator_string, 
    Correlator_list& correlator_list){

  for(auto str : correlator_string){
  
    std::vector<std::string> correlator_tokens;
    boost::split(correlator_tokens, str, boost::is_any_of(":"));
  
    std::string type;
    std::vector<int> quark_number;
    std::vector<int> operator_number;
    std::string GEVP;
    std::vector<int> tot_mom;
    for (auto corr_t : correlator_tokens){
      // getting the type name
      if (corr_t.compare(0,1,"C") == 0)
        type = corr_t;
      // getting quark numbers
      else if (corr_t.compare(0,1,"Q") == 0) 
        quark_number.push_back(boost::lexical_cast<int>(corr_t.erase(0,1)));
      // getting operator numbers
      else if (corr_t.compare(0,2,"Op") == 0)
        operator_number.push_back(boost::lexical_cast<int>(corr_t.erase(0,2)));
      // getting the GEVP type
      else if (corr_t.compare(0,1,"G") == 0)
        GEVP = corr_t;
      // getting total momenta for moving frames
      else if (corr_t.compare(0,1,"P") == 0) {
        corr_t.erase(0,1);
        std::vector<std::string> tokens;
        boost::split(tokens, corr_t, boost::is_any_of(","));
        for(auto t : tokens)
          tot_mom.push_back(boost::lexical_cast<int>(t));
      }
      // catching wrong entries
      else {
        std::cout << "there is something wrong with the correlators" 
                  << std::endl;
        exit(0);
      }
    }
    correlator_list.push_back(Correlators
                          (type, quark_number, operator_number, GEVP, tot_mom));
  }
  // TODO: write check for correctness of input data
}

// *****************************************************************************
// simplifies and cleans read_parameters function
void config_input_data_handling (const int start_config,
    const int end_config, const int delta_config) {

  try{
    if(start_config < 0){
      std::cout << "\ninput file error:\n" << "\toption \"start config\""
          << " is mandatory and its value must be an integer greater or equal 0!"
          << "\n\n";
      exit(0);
    }
    else if(end_config < 1 || end_config < start_config){
      std::cout << "\ninput file error:\n" << "\toption \"end_config\""
          << " is mandatory, its value must be an integer greater than 0,"
          << " and it must be larger than start config!" << "\n\n";
      exit(0);
    }
    else if(delta_config < 1){
      std::cout << "\ninput file error:\n" << "\toption \"delta_config\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tprocessing configurations " << start_config << " to "
        << end_config << " in steps of " << delta_config << "\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

} // end of unnamed namespace

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void GlobalData::input_handling(
    const std::vector<std::string>& quark_configs,
    const std::vector<std::string>& operator_list_configs,
    const std::vector<std::string>& correlator_list_configs) {

  lattice_input_data_handling(path_output, name_lattice, 
                                                path_config, Lt, Lx, Ly, Lz);
  eigenvec_perambulator_input_data_handling(
      number_of_eigen_vec, path_eigenvectors, name_eigenvectors, 
      path_perambulators, name_perambulators);
  quark_input_data_handling(quark_configs, quarks);
  operator_input_data_handling(operator_list_configs, operator_list);
  correlator_input_data_handling(correlator_list_configs, correlator_list);
  config_input_data_handling(start_config, end_config, 
                                               delta_config);
  // TODO: Are these still needed anywhere?
  // computing some global variables depending on the input values
  dim_row = Lx * Ly * Lz * 3;

  //needed for config_utils.h
  //4 is number of directions, 3 number of colors and 2 factor
  //for memory requirement of complex numbers
  V_TS = dim_row * 4 * 3 * 2;
  V_for_lime = V_TS * Lt;

}

