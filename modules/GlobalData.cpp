/*
 * GlobalData.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#include "GlobalData.h"

namespace po = boost::program_options;

GlobalData* GlobalData::instance_ = 0;

GlobalData* GlobalData::Instance () {

  if(instance_ == 0) instance_ = new GlobalData;

  return instance_;
}
// *****************************************************************************
// A helper function to simplify the main part.
template<class T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}
// *****************************************************************************
/// @brief Convenience function for when a 'store_to' value is being provided
///        to typed_value.
///
/// @param store_to The variable that will hold the parsed value upon notify.
///
/// @return Pointer to a type_value.
template<typename T>
boost::program_options::typed_value<T>* make_value (T* store_to) {
  return boost::program_options::value<T>(store_to);
}
// *****************************************************************************
/// @brief Stream insertion operator for slave.
///
/// @param stream The stream into which quark is being inserted.
/// @param q The quark object.
///
/// @return Reference to the ostream.
static std::ostream& operator<< (std::ostream& stream, const quark& quark) {
  return stream << "\tQUARK type: ****  " << quark.type
      << "  ****\n\t number of random vectors: " << quark.number_of_rnd_vec
      << "\n\t dilution scheme in time: " << quark.dilution_T
      << quark.number_of_dilution_T << "\n\t dilution scheme in ev space: "
      << quark.dilution_E << quark.number_of_dilution_E
      << "\n\t dilution scheme in Dirac space: " << quark.dilution_D
      << quark.number_of_dilution_D << "\n";
}
// *****************************************************************************
/// @brief Makes a quark object from a string
quark make_quark (const std::string& quark_string) {
  // Tokenize the string on the ":" delimiter.
  std::vector<std::string> tokens;
  boost::split(tokens, quark_string, boost::is_any_of(":"));

  // If the split did not result in exactly 8 tokens, then the value
  // is formatted wrong.
  if(8 != tokens.size()){
    using boost::program_options::validation_error;
    throw validation_error(validation_error::invalid_option_value,
        "quarks.quark", quark_string);
  }

  // Create a quark from the token values.
  return quark(tokens[0], boost::lexical_cast<int>(tokens[1]), tokens[2],
      boost::lexical_cast<int>(tokens[3]), tokens[4],
      boost::lexical_cast<int>(tokens[5]), tokens[6],
      boost::lexical_cast<int>(tokens[7]));
}

////////////////////////////////////////////////////////////////////////////////
//TODO: clean that up ////////////////////////////////////////////////////////// 
//init operator ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void GlobalData::init_from_infile() {

  //TODO: think about the numbers of momenta, displacements and gamma
  const size_t nb_mom = get_number_of_momenta(); //momentum_squared.size();
  const size_t nb_mom_sq = number_of_momentum_squared;
  //TODO: include displacement into dg (displacementgamma) multiindex
  const size_t nb_dg = number_of_displ_gamma;

  // nb_op - number of combinations of three-momenta and gamma structures
  // op    - vector of all three-momenta, three-displacements and gamma 
  //         structure combinations
  const size_t nb_op = number_of_operators;
  op_Corr.resize(nb_op);
  set_Corr();

  // nb_op_C2 - number of combinations of absolute values squared of momenta
  //            and gamma-displacement combinations for 2pt-fct
  // op_C2    - vector of all combinations for 2pt-fct and vector of 
  //            op-index-pairs with all corresponding three-vectors and gammas
  const size_t nb_op_C2 = nb_mom_sq * nb_dg * nb_dg;
  op_C2.resize(nb_op_C2);
  set_C2();

  const size_t nb_op_C4 = nb_mom_sq * nb_dg * nb_dg;
  op_C4.resize(nb_op_C4);
  set_C4();

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// function to initialize the vector of all necessary operators
void GlobalData::set_Corr(){

  const int max_mom_squared = number_of_max_mom;
  const size_t nb_dg = dg.size();

  size_t i = 0;

  // all three-momenta
  for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
    for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
      for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
        if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
          continue;
        }
        // all gamma structures
        for(size_t gam = 0; gam < nb_dg; gam++){

          op_Corr[i].p[0] = ipx;
          op_Corr[i].p[1] = ipy;
          op_Corr[i].p[2] = ipz;
          op_Corr[i].dis[0] = 0;
          op_Corr[i].dis[1] = 0;
          op_Corr[i].dis[2] = 0;
          op_Corr[i].gamma[0] = dg[gam];
          op_Corr[i].gamma[1] = 4;
          op_Corr[i].gamma[2] = 4;
          op_Corr[i].gamma[3] = 4;
  
          op_Corr[i].id = i;
          i++;
        }
      }
    }
  }

  if(i != op_Corr.size()){
    std::cout << "Error in LapH::set_Cor(): nb_op not equal to allocated "
                 "number of operators" << std::endl;
    exit(0);
  } 

}

// function to obtain the index combinations in op for the 2pt-fct
void GlobalData::set_C2(){

  const size_t nb_mom_sq = number_of_momentum_squared;
  const size_t nb_dg = dg.size();

  size_t nb_op = op_Corr.size();
  size_t j = 0;

  // run over all momenta squared (back-to-back hardcoded) and gamma 
  // combinations
  for(size_t p_sq = 0; p_sq < nb_mom_sq; p_sq++){
    for(size_t so = 0; so < nb_dg; so++){
      for(size_t si = 0; si < nb_dg; si++){

        // index for access of element
        size_t i = p_sq * nb_dg*nb_dg + so * nb_dg + si;

        // save p^2 and gamma structure at source and sink
        op_C2[i].p_sq = p_sq;
        op_C2[i].dg_so = so;
        op_C2[i].dg_si = si;

        // loop over op and set index pairs
        for(auto& el : op_Corr)
          if((el.p[0]*el.p[0] + el.p[1]*el.p[1] + el.p[2]*el.p[2]) == p_sq){ 
            if(el.gamma[0] == dg[so]){
              size_t id1 = el.id;
              // thats the generalized version of nb_mom - p - 1 including 
              // a faster running gamma structure
              size_t id2 = nb_op - nb_dg * (id1/nb_dg + 1) + si;
              // warning because array has no list-of-arrays constructor but 
              // works. Can change this to pair structure.
              op_C2[i].index.emplace_back(std::pair<size_t, size_t>(id1, id2));
            }
          }

        j++;
        
      }
    }
  }

  if(j != op_C2.size()){
    std::cout << "Error in LapH::set_C2(): nb_op not equal to allocated "
                 "number of operators" << std::endl;
    exit(0);
  } 

}

// function to obtain the index combinations in op for the 2pt-fct
void GlobalData::set_C4(){

  const size_t nb_mom_sq = number_of_momentum_squared;
  const size_t nb_dg = dg.size();

  size_t nb_op = op_Corr.size();
  size_t j = 0;

  // run over all momenta squared (back-to-back hardcoded) and gamma 
  // combinations
  for(size_t p_sq = 0; p_sq < nb_mom_sq; p_sq++){
    //only diagonal elements
    size_t p_sq_so = p_sq;
    size_t p_sq_si = p_sq;
    for(size_t so = 0; so < nb_dg; so++){
    for(size_t si = 0; si < nb_dg; si++){

      // index for access of element
      size_t i = p_sq*nb_dg*nb_dg + so*nb_dg + si;

      // save p^2 and gamma structure at source and sink
      op_C4[i].p_sq_so = p_sq_so;
      op_C4[i].p_sq_si = p_sq_si;
      op_C4[i].dg_so = so;
      op_C4[i].dg_si = si;

      // loop over op and set index pairs
      for(auto& el_so : op_Corr)
        if((el_so.p[0]*el_so.p[0] + el_so.p[1]*el_so.p[1] + 
            el_so.p[2]*el_so.p[2]) == p_sq_so){ 
        if(el_so.gamma[0] == dg[so]){
          size_t id1 = el_so.id;
          // thats the generalized version of nb_mom - p - 1 including 
          // a faster running gamma structure
          size_t id2 = nb_op - nb_dg * (id1/nb_dg + 1) + so;

          for(auto& el_si : op_Corr)
            if((el_si.p[0]*el_si.p[0] + el_si.p[1]*el_si.p[1] + 
                el_si.p[2]*el_si.p[2]) == p_sq_si){ 
            if(el_si.gamma[0] == dg[si]){
              size_t id3 = el_si.id;
              size_t id4 = nb_op - nb_dg * (id3/nb_dg + 1) + si;

              op_C4[i].index.emplace_back(
                  std::array<size_t, 4>{{id1, id2, id3, id4}});
            }}//loops over sink
        }}//loops over source

      j++;
      
    }}//loops over displ-gamma
  }//loop over mom_sq

  if(j != op_C4.size()){
    std::cout << "Error in LapH::set_C4(): nb_op not equal to allocated "
                 "number of operators" << std::endl;
    exit(0);
  } 

}
// *****************************************************************************
// simplifies and cleans read_parameters function
static void lattice_input_data_handling (const std::string path_output,
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
static void eigenvec_perambulator_input_data_handling (
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
// *****************************************************************************
// simplifies and cleans read_parameters function
static void momentum_input_data_handling (const int number_of_max_mom,
    const int max_mom_in_one_dir, std::vector<int>* mom_squared) {

  try{
    if(number_of_max_mom < 0){
      std::cout << "\ninput file error:\n" << "\toption \"number_of_max_mom\""
          << " is mandatory and its value must be an integer greater or equal 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tabsolute value squared of max momentum .... "
        << number_of_max_mom << "\n";
    if(max_mom_in_one_dir < 0){
      std::cout << "\ninput file error:\n" << "\toption \"max_mom_in_one_dir\""
          << " is mandatory and its value must be an integer greater or equal 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tmaximal momentum in one direction ........ "
        << max_mom_in_one_dir << "\n";

    int max_mom_squared = number_of_max_mom;
    // generate all used momenta
    for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
      for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
        for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
          if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
            continue;
          }
          mom_squared->push_back(ipx * ipx + ipy * ipy + ipz * ipz);
          std::cout << "\tmomentum p = " << mom_squared->size() - 1
            << " corresponds to ............ (" 
            << ipx << ", " << ipy << ", " << ipz << ")\n" << std::endl;
        }
      }
    }

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
static void dirac_input_data_handling (const int dirac_min,
    const int dirac_max) {

  try{
    if(dirac_min < 0 || dirac_min > 15){
      std::cout << "\ninput file error:\n" << "\toption \"dirac_min\""
          << " is mandatory and its value must be an integer greater or equal 0 and smaller 16!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tlowest Dirac index used .................. "
        << dirac_min << "\n";
    if(dirac_max < 0 || dirac_max > 15){
      std::cout << "\ninput file error:\n" << "\toption \"dirac_max\""
          << " is mandatory and its value must be an integer greater or equal 0 and smaller 16!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\thighest Dirac index used ................. "
        << dirac_max << "\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
static void displacement_input_data_handling (const int displ_min,
    const int displ_max) {

  try{
    if(displ_min < 0 || displ_min > 3){
      std::cout << "\ninput file error:\n" << "\toption \"displ_min\""
          << " is mandatory and its value must be an integer greater or equal 0 and smaller 4!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tmimimal displacement used ................ "
        << displ_min << "\n";
    if(displ_max < 0 || displ_max > 3){
      std::cout << "\ninput file error:\n" << "\toption \"displ_max\""
          << " is mandatory and its value must be an integer greater or equal 0 and smaller 4!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tmaximal displacement used ................ "
        << displ_max << "\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
static void config_input_data_handling (const int start_config,
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
// *****************************************************************************
// simplifies and cleans read_parameters function
static void quark_check (quark quarks) {

  try{
    if(quarks.type != "u" && quarks.type != "d" && quarks.type != "s"
        && quarks.type != "c"){
      std::cout << "quarks.quark.type must be u, d, s or c" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_rnd_vec < 1){
      std::cout << "quarks.quark.number_of_rnd_vec must be greater than 0"
          << std::endl;
      exit(0);
    }
    else if(quarks.dilution_T != "TI" && quarks.dilution_T != "TB"){
      std::cout << "quarks.quark.dilutione_T must be TI or TB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_T < 1){
      std::cout << "quarks.quark.number_of_dilution_T must be greater than 0 "
          "and smaller than the temporal extend" << std::endl;
      exit(0);
    }
    else if(quarks.dilution_E != "EI" && quarks.dilution_E != "EB"){
      std::cout << "quarks.quark.dilutione_E must be EI or EB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_E < 1){
      std::cout << "quarks.quark.number_of_dilution_E must be greater than 0 "
          "and smaller than number of eigen vectors" << std::endl;
      exit(0);
    }
    else if(quarks.dilution_D != "DI" && quarks.dilution_D != "DI"){
      std::cout << "quarks.quark.dilutione_D must be DI or DB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_D < 1 || quarks.number_of_dilution_D > 4){
      std::cout << "quarks.quark.number_of_dilution_D must be greater than 0 "
          "and smaller than 5" << std::endl;
      exit(0);
    }
    else std::cout << quarks << std::endl;
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }

}
void GlobalData::quark_input_data_handling (
    const std::vector<std::string> quark_configs) {
  try{
    // Transform each configured quark into a quark via make_quark, inserting each
    // object into the quark vector.
    std::transform(quark_configs.begin(), quark_configs.end(),
        std::back_inserter(quarks), make_quark);
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
void GlobalData::read_parameters (int ac, char* av[]) {

  try{
    std::string input_file;
    std::string output_file;
    // Variables that will store parsed values for quarks.
    std::vector<std::string> quark_configs;

    // Declare a group of options that will be allowed only on command line
    po::options_description generic("Command line options");
    generic.add_options()("help,h", "produce help message")("version,v",
        "print version string")("verbose",
        "does additional tests and prints more details")("input,i",
        po::value<std::string>(&input_file)->default_value("LapHs.in"),
        "name of input file.")("output,o",
        po::value<std::string>(&output_file)->default_value("LapHs.out"),
        "name of output file.");

    // Declare a group of options that will be
    // allowed both on command line and in input file
    po::options_description config("Input file options");
    // lattice options
    config.add_options()("output_path",
        po::value<std::string>(&path_output)->default_value("../../contractions"),
        "path for output")
        ("config_path",
        po::value<std::string>(&path_config)->default_value("../../configs"),
        "path for configurations")
        ("lattice", po::value<std::string>(&name_lattice)->
        default_value("lattice"),"Codename of the lattice")("Lt", 
        po::value<int>(&Lt)->default_value(0),
        "Lt: temporal lattice extend")("Lx",
        po::value<int>(&Lx)->default_value(0),
        "Lx: lattice extend in x direction")("Ly",
        po::value<int>(&Ly)->default_value(0),
        "Ly: lattice extend in y direction")("Lz",
        po::value<int>(&Lz)->default_value(0),
        "Lz: lattice extend in z direction");
    // eigenvector options
    config.add_options()("number_of_eigen_vec",
        po::value<int>(&number_of_eigen_vec)->default_value(0),
        "Number of eigen vectors")("path_eigenvectors",
        po::value<std::string>(&path_eigenvectors)->default_value("."),
        "directory of eigenvectors")("name_eigenvectors",
        po::value<std::string>(&name_eigenvectors)->default_value(
            "eigenvector"),
        "name of eigenvectors\nThe full name is internally created to:\n"
            "\"name_of_eigenvectors.eigenvector\n."
            "time slice.configuration\"");
    // perambulator options
    config.add_options()("path_perambulators",
        po::value<std::string>(&path_perambulators)->default_value("."),
        "directory of perambulators")("name_perambulators",
        po::value<std::string>(&name_perambulators)->default_value(
            "perambulator"),
        "name of perambulators\nThe full name is internally created to:\n"
            "\"rather long\n."
            "t_sink.configuration\"");
    // quark options
    config.add_options()("quarks.quark", make_value(&quark_configs),
        "quark input, must be of type:\n"
            "quark = \n type:number of rnd. vec.:\n"
            " dil type time:number of dil time:\n"
            " dil type ev:number of dil ev:\n"
            " dil type Dirac:number of dil Dirac");
    // momentum options
    config.add_options()("number_of_max_mom",
        po::value<int>(&number_of_max_mom)->default_value(-1),
        "Maximum momentum squared")("max_mom_in_one_dir",
        po::value<int>(&max_mom_in_one_dir)->default_value(-1),
        "Maximum momentum in one direction");
    // dirac options
    config.add_options()("dirac_min",
        po::value<int>(&dirac_min)->default_value(-1),
        "dirac_min")("dirac_max",
        po::value<int>(&dirac_max)->default_value(-1),
        "dirac_max");
    // displacement options
    config.add_options()("displ_min",
        po::value<int>(&displ_min)->default_value(-1),
        "displ_min")("displ_max",
        po::value<int>(&displ_max)->default_value(-1),
        "displ_max");
    // configuration options
    config.add_options()("start_config",
        po::value<int>(&start_config)->default_value(-1), "First configuration")(
        "end_config", po::value<int>(&end_config)->default_value(0),
        "Last configuration")("delta_config",
        po::value<int>(&delta_config)->default_value(0),
        "Stepsize between two configurations");

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description input_file_options;
    input_file_options.add(config);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);
    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    po::store(
        po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(),
        vm);
    po::notify(vm);

    // command line options ****************************************************
    if(vm.count("help")){
      std::cout << visible << "\n";
      exit(0);
    }
    if(vm.count("verbose")){
      verbose = 1;
    }
    else verbose = 0;
    if(vm.count("version")){
      std::cout << "stochastic LapH code, version under construction \n";
      exit(0);
    }
    std::ifstream ifs(input_file.c_str());
    if(!ifs){
      std::cout << "CANNOT open input file: " << input_file << "\n";
      exit(0);
    }
    else{
      po::store(parse_config_file(ifs, input_file_options), vm);
      po::notify(vm);
    }
    ifs.close();

    // input file options ******************************************************
    //
    lattice_input_data_handling(path_output, name_lattice, path_config, Lt, Lx, Ly, Lz);
    //
    eigenvec_perambulator_input_data_handling(number_of_eigen_vec,
        path_eigenvectors, name_eigenvectors, path_perambulators,
        name_perambulators);
    //
    quark_input_data_handling(quark_configs);
    //
    momentum_input_data_handling(number_of_max_mom, max_mom_in_one_dir, &momentum_squared);
    //
    dirac_input_data_handling(dirac_min, dirac_max);
    //
    displacement_input_data_handling(displ_min, displ_max);
    //
    config_input_data_handling(start_config, end_config, delta_config);

    // computing some global variables depending on the input values ***********
    dim_row = Lx * Ly * Lz * 3;

    //needed for config_utils.h
    //4 is number of directions, 3 number of colors and 2 factor
    //for memory requirement of complex numbers
    V_TS = dim_row * 4 * 3 * 2;
    V_for_lime = V_TS * Lt;

    //dirac structure hard coded
    dg.push_back(5);

    //displacement not supported yet
    number_of_displ_gamma = dg.size();
    number_of_operators = momentum_squared.size() * number_of_displ_gamma;
    number_of_momentum_squared = number_of_max_mom + 1;

    init_from_infile();

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
