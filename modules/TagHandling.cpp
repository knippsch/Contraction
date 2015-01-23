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

// create array for gamma from stirng
static std::array<int, 4> create_gamarray_from_string(std::string in) {
  std::array<int, 4> out;
  std::vector<std::string> tokens;
  // erasing the brakets at the beginning and the end 
  boost::split(tokens, in, boost::is_any_of("g"));
  tokens.erase(tokens.begin());
  //for (auto& el : tokens) std::cout << el << std::endl;
  for (size_t ind_out = 0; ind_out < tokens.size(); ++ind_out){
    out[ind_out] = boost::lexical_cast<int>(tokens[ind_out]); 
  }
  // fill 0 in for the rest
  for (size_t rest = tokens.size(); rest < 4; ++rest ) out[rest] = 0;
  return out;
}

//  turn string into tag for use in CorrelatorIo.cpp
void string_to_tag(const std::string& in, Tag& out){
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

  boost::split(splits, in, boost::is_any_of(":"));
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
    std::string g_str;
    for(auto& qu : pos){
      //sort strings according to type
      if(qu.front() == 'p') p_str.push_back(qu);
      else if(qu.front() == 'g') g_str += qu;
      else if(qu.front() == 'd') d_str.push_back(qu);
    }

      for(auto& el : p_str) out.mom.push_back(create_3darray_from_string(el));
      for(auto& el : d_str) out.dis.push_back(create_3darray_from_string(el));
      out.gam.push_back(create_gamarray_from_string(g_str));
  }
  out.mom_cm = square_comp(out.mom[0], out.mom[1]);
  for(auto& el : quarks) out.q_cont += el;
   
}

void tag_to_string(const Tag& in, std::string& out){
  size_t no_part = in.q_cont.length();
  std::cout << in.q_cont << std::endl;
  for(size_t i = 0; i < no_part; ++i){
    //get string of quantum each quantum number
    out.push_back(in.q_cont[i]);
    out += ":";
    //out += (std::c_str(in.q_cont[i]) + ":");
    // gammas of each point go first
    out += "g";
    std::vector<std::string> tmp;
    for(auto& el : in.gam[i]) tmp.push_back(std::to_string(el));
    std::string joined = boost::algorithm::join(tmp,".g");
    out += (joined + ".");
    tmp.clear();

    //then displacement
    for(auto& el : in.dis[i]) tmp.push_back(std::to_string(el));
    joined = boost::algorithm::join(tmp,",");
    out += "d(" + joined + ").";
    tmp.clear();

    //at last momentum
    for(auto& el : in.mom[i]) tmp.push_back(std::to_string(el));
    joined = boost::algorithm::join(tmp,",");
    out += "p(" + joined +"):";
  }
}

// print tag
void print_tag(const Tag& query){
  std::cout << "Quark content: " << query.q_cont << std::endl;
  std::cout << "Center of mass momentum is: " << query.mom_cm  << std::endl;
  std::cout << "momenta:" << std::endl;
  for (auto& el : query.mom ) {
    for (auto& comp : el) std::cout << comp << " ";
    std::cout << " " << std::endl;
  }
   std::cout << " \n displacements:" << std::endl;
  for (auto& el : query.dis ) {
    for (auto& comp : el) std::cout << comp << " ";
    std::cout << " " << std::endl;
  }
   std::cout << " \n gammas:" << std::endl;
  for (auto& el : query.gam ) {
    for (auto& comp : el) std::cout << comp << " ";
    std::cout << " " << std::endl;
  }
}

// Compare two tags of correlation functions
bool compare_tags(const Tag& tag1, const Tag& tag2){
  bool flag = true;
  if (tag1.q_cont != tag2.q_cont) flag = false;
  if (tag1.mom_cm != tag2.mom_cm) flag = false;
  for(size_t p = 0; p < tag1.mom.size(); ++p){
    if (tag1.mom[p] != tag2.mom[p]) flag = false;
    if (tag1.dis[p] != tag2.dis[p]) flag = false;
    if (tag1.gam[p] != tag2.gam[p]) flag = false;
  }
  return flag;
}


void zero_tag(const size_t pts, Tag& def){
  def.mom_cm = 0;
  std::array<int, 3 > zero_3d {{0,0,0}};
  std::array<int, 4 > zero_4d {{0,0,0,0}};
  for (size_t ind = 0; ind < pts; ++ind) {
    def.mom.push_back(zero_3d);
    def.dis.push_back(zero_3d);
    def.gam.push_back(zero_4d);
  } 
}

void zero_vec_tag(const size_t pts, std::vector<Tag> attr){
  for (auto& el : attr) zero_tag(pts, el);
}













