#include "global_data.h"
#include "global_data_utils.h"

namespace {

// *****************************************************************************
// A helper function to simplify the main part.
template<class T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}
// *****************************************************************************
/// @brief Stream insertion operator for slave.
///
/// @param stream The stream into which quark is being inserted.
/// @param q The quark object.
///
/// @return Reference to the ostream.
std::ostream& operator<< (std::ostream& stream, const quark& quark) {
  return stream << "\tQUARK type: ****  " << quark.type
      << "  ****\n\t number of random vectors: " << quark.number_of_rnd_vec
      << "\n\t dilution scheme in time: " << quark.dilution_T
      << quark.number_of_dilution_T << "\n\t dilution scheme in ev space: "
      << quark.dilution_E << quark.number_of_dilution_E
      << "\n\t dilution scheme in Dirac space: " << quark.dilution_D
      << quark.number_of_dilution_D << "\n";
}


// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
std::array<int, 3> create_3darray_from_string(std::string in) { 

  std::array<int, 3> out;
  std::vector<std::string> tokens;
  // erasing the brakets at the beginning and the end
  in.erase(0,2);
  in.erase(in.end()-1);

  boost::split(tokens, in, boost::is_any_of(","));

  return {{boost::lexical_cast<int>(tokens[0]),
          boost::lexical_cast<int>(tokens[1]),
          boost::lexical_cast<int>(tokens[2]) }};

}
// *****************************************************************************
void create_all_momentum_combinations(const int p, 
                                        std::vector<std::array<int, 3> >& out) {
  // creating all momentum combinations possible and needed
  int max_p = p;
  std::vector<std::array<int, 3> > all_p;
  for(int p1 = -max_p; p1 < max_p+1; p1++)
    for(int p2 = -max_p; p2 < max_p+1; p2++)
      for(int p3 = -max_p; p3 < max_p+1; p3++)
        all_p.push_back({{p1, p2, p3}});
  // copying wanted combinations into out array
  for(const auto& all : all_p)
    if(p == all[0]*all[0] + all[1]*all[1] + all[2]*all[2])
      out.push_back(all);

}
// *****************************************************************************
void create_mom_array_from_string(std::string in, 
                                        std::vector<std::vector<std::array
                                                   <int, 3> > >& out) {
  // erase the p (first entry)
  in.erase(0,1);
  std::vector<std::string> tokens;
  boost::split(tokens, in, boost::is_any_of(","));
  int p;
  size_t counter = 0;
  out.resize(tokens.size());
  for(const auto& t : tokens){
    p = boost::lexical_cast<int>(t);
    create_all_momentum_combinations(p, out[counter]);
    counter++;
  }

}

} // end of unnamed namespace

namespace global_data_utils {

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

// *****************************************************************************
// simplifies and cleans read_parameters function
void quark_check (quark quarks) {

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

// *****************************************************************************
/// @brief Makes an operator list object from a string
Operator_list make_operator_list(const std::string& operator_string) {

  Operator_list op_list; // return object

  // Two steps are necessary: 1. Getting all operators in one list which are 
  //                             separated by ";"
  //                          2. Separating the individual operators into its
  //                             smaller bits, which are separated by "."
  // Tokenize the string on the ";" delimiter -> Individual operators
  std::vector<std::string> operator_tokens;
  boost::split(operator_tokens, operator_string, boost::is_any_of(":"));

  // running over opeator tokens and split them further (Step 2):
  for (const auto& op_t : operator_tokens){
    std::vector<std::string> tokens;
    boost::split(tokens, op_t, boost::is_any_of("."));
    std::vector<int> gammas;
    std::array<int, 3> dil_vec;
    std::vector<std::vector<std::array<int, 3> > >mom_vec;
    for (auto str : tokens){
      // getting the gamma structure
      if(str.compare(0,1,"g") == 0)
        gammas.push_back(boost::lexical_cast<int>(str.erase(0,1)));
      // getting the displacement indices
      else if (str.compare(0,1,"d") == 0) {
        if(str.compare(1,1,"0") == 0)
          dil_vec = {{0, 0, 0}};
        else if (str.compare(1,1,"(") == 0)
          dil_vec = create_3darray_from_string(str);
        else {
         std::cout << "Something wrong with the displacement in the operator" \
                      " definition" << std::endl;
         exit(0);
        }
      }
      // getting the momenta
      else if (str.compare(0,1,"p") == 0) {
        if(str.compare(1,1,"(") == 0){
          mom_vec.resize(1);
          mom_vec[0].push_back(create_3darray_from_string(str));
        }
        else 
          create_mom_array_from_string(str, mom_vec);
      }
      // catching wrong entries
      else {
        std::cout << "there is something wrong with the operators" << std::endl;
        exit(0);
      }
    }
    op_list.push_back(Operators(gammas, dil_vec, mom_vec));
  }
  return op_list;
}

} // end of namespace global_data_utils
