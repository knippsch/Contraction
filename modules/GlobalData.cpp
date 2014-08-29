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
/// @brief Makes a slave given an address and port.
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

// *****************************************************************************
// simplifies and cleans read_parameters function
static void lattice_input_data_handling (const std::string path_output,
    const std::string name_lattice, int Lt, int Lx, int Ly, int Lz) {
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
		std::cout << "\tResults will be saved to path:\n\t\t"
				<< path_output << "/" << name_lattice << "/\n";
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
		const int max_mom_in_one_dir) {

	try{
		if(number_of_max_mom < 0){
			std::cout << "\ninput file error:\n" << "\toption \"number_of_max_mom\""
					<< " is mandatory and its value must be an integer greater or equal 0!"
					<< "\n\n";
			exit(0);
		}
		else std::cout << "\ttotal number of momenta .................. "
				<< number_of_max_mom << "\n";
		if(max_mom_in_one_dir < 0){
			std::cout << "\ninput file error:\n" << "\toption \"max_mom_in_one_dir\""
					<< " is mandatory and its value must be an integer greater or equal 0!"
					<< "\n\n";
			exit(0);
		}
		else std::cout << "\tmaximal momentum in one direction ........ "
				<< max_mom_in_one_dir << "\n\n";
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
		else std::cout << "\tdirac_min = "
				<< dirac_min << "\n";
		if(dirac_max < 0 || dirac_max > 15){
			std::cout << "\ninput file error:\n" << "\toption \"dirac_max\""
					<< " is mandatory and its value must be an integer greater or equal 0 and smaller 16!"
					<< "\n\n";
			exit(0);
		}
		else std::cout << "\tdirac_max = "
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
		else std::cout << "\tdispl_min = "
				<< displ_min << "\n";
		if(displ_max < 0 || displ_max > 3){
			std::cout << "\ninput file error:\n" << "\toption \"displ_max\""
					<< " is mandatory and its value must be an integer greater or equal 0 and smaller 4!"
					<< "\n\n";
			exit(0);
		}
		else std::cout << "\tdispl_max = "
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
		else if(end_config < 0 || end_config < start_config){
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
		lattice_input_data_handling(path_output, name_lattice, Lt, Lx, Ly, Lz);
		//
		eigenvec_perambulator_input_data_handling(number_of_eigen_vec,
				path_eigenvectors, name_eigenvectors, path_perambulators,
				name_perambulators);
		//
		quark_input_data_handling(quark_configs);
		//
		momentum_input_data_handling(number_of_max_mom, max_mom_in_one_dir);
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
	}
	catch(std::exception& e){
		std::cout << e.what() << "\n";
		exit(0);
	}
}
