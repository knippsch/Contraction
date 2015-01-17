// This file contains the whole functionality for handling the tags of the
// correlation function outputs
// Idea is to use a structure similar for the infile for the handling of the
// Metadata used for the contraction, the strings should be written into the lime
// record

#include "TagHandling.h"

// compose string to look for in the outfile from infile
void compose_string(const char* filename, std::string& tag){
  // substrings needed for search string
  std::string config_line;
  // strings storing identifier and value
  std::string entry;
  std::string option;
  // sort types into list
  std::string corr_type;
  std::string misc;
  std::vector<std::string> quarks;
  std::vector<std::string> qu_numbers;
  // open file containing quantum numbers
  std::fstream quanta;
  quanta.open(filename, std::fstream::in);
  // get lines
  while (!std::getline(quanta, config_line).eof()){
    if(config_line[0] != '#') {
      //find first occurence of '=' in each line and divide strings
      unsigned pos = config_line.find_first_of("=");
      entry = config_line.substr(0,(pos-1));
      option = config_line.substr( (pos+2), std::string::npos);
      if(entry=="function") corr_type = option;
      if(entry=="operator_list") qu_numbers.push_back(option);
      if(entry=="quark") quarks.push_back(option);
      if(entry=="misc") misc = option;
    }
  }
  // all is read in, file can be closed
  quanta.close();
  // put everything together as search string
  // Check if corr_type complies with operator number
  size_t num_quarks = stoi(corr_type.substr(1,2));
  if(quarks.size()!=num_quarks){
    std::cout << "Wrong number of quarks " << corr_type << std::endl;
    std::cout << quarks.size()<< " given, " << num_quarks << " expected" <<
      std::endl;
  }
  else{
    tag+=corr_type;
    tag+=":";
    for(size_t num_op = 0; num_op < quarks.size(); ++num_op){
      tag+=quarks[num_op];
      tag+=":";
      tag+=qu_numbers[num_op];
      tag+=":";
    }
    tag+=misc;
  }
}
// create array from stirng
static std::array<int, 3> create_3darray_from_string(std::string in) {
  std::array<int, 3> out;
  std::vector<std::string> tokens;
  // erasing the brakets at the beginning and the end
  if (in[1] == '('){
    in.erase(0,2);
    in.erase(in.end()-1);
  }
  else{
    in.erase(0,1);
  }
  boost::split(tokens, in, boost::is_any_of(","));
  if(tokens.size() == 3){
  return {{boost::lexical_cast<int>(tokens[0]),
    boost::lexical_cast<int>(tokens[1]),
    boost::lexical_cast<int>(tokens[2]) }};
  }
  else{
    return{{0, 0, 0}};
  }
}

//  turn string into tag for use in CorrelatorIo.cpp
void string_to_tag(const std::string& search, Tag sign){
  std::vector<std::string> splits;

  // Quantum numbers
  std::string particle;
  std::string misc;
  std::vector<std::string> quarks;
  std::vector<std::string> operat;
  std::vector<std::vector<std::string > > quanta;

  // quantum numbers as arrays
  std::vector<std::string > mom;
  std::vector<std::string > dis;
  std::vector<std::string > gam;

  // quantum numbers as vectors and ints
  std::vector<std::array<int, 3 > > mom_i;
  std::vector<std::array<int, 3 > > dis_i;
  std::vector<int> gam_i;

  boost::split(splits, search, boost::is_any_of(":"));
  particle = splits.front();
  misc = splits.back();
  splits.erase(splits.begin());
  splits.pop_back();
  // sort into quarks and operators
  for(auto& c : splits) {
    if(c.size() == 1) quarks.push_back(c);
    else operat.push_back(c);
  }
  // each particle gets an own set of vectors for its quantum numbers
  for(auto& op : operat){
    std::vector<std::string > particle;
    split(particle, op, boost::is_any_of("."));
    quanta.push_back(particle);
  } 
  // loop over all particles
  for(auto& pos : quanta){
    // each particle has 3 types of attributes
    std::vector<std::string > p_str;
    std::vector<std::string > d_str;
    std::vector<std::string > g_str;
    // same in an integer version
    std::vector<std::array<int, 3 > > p_int;
    std::vector<std::array<int, 3 > > d_int;
    std::vector<int> g_int;
    for(auto& qu : pos){
    //sort strings according to type
        //std::cout << qu << std::endl;
        if(qu.front() == 'p') p_str.push_back(qu);
        else if(qu.front() == 'g') g_str.push_back(qu);
        else if(qu.front() == 'd') d_str.push_back(qu);
    }
    for(auto& el : p_str) sign.mom.push_back(create_3darray_from_string(el));
    for(auto& el : d_str) sign.dis.push_back(create_3darray_from_string(el));
    for(auto& el : g_str) sign.gam.push_back((boost::lexical_cast<int>(el.substr(1))));
 } 
}

















