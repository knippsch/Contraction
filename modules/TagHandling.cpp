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
      if(entry=="miscellanous") misc = option;
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

//  turn string into tag for use in CorrelatorIo.cpp
void string_to_tag(const std::string& search, Tag sign){

  // bookkeeping variables
  std::string tmp;
  std::vector<size_t> pos_col;

  // Quantum numbers
  std::string particle;
  std::vector<std::string> quarks;
  std::vector<std::string> operat;

  // quantum numbers as arrays
  std::vector<std::array<int, 3 > > mom;
  std::vector<std::array<int, 3 > > dis;
  std::vector<size_t > gam;
  
  // get all positions of delimiters
  for(size_t pos = 0; pos < search.length(); ++pos){
    char c = search[pos];
    if(c == ':') pos_col.push_back(pos);
  }
  for(auto& c: pos_col) std::cout << c << std::endl;
  // first entry always is correlator type
  particle = search.substr(0, pos_col[0]);
  std::cout << "particle type is: " << particle << std::endl;
  // sort quarks and operators
  //  loop has to end earlier, want to compare distances
  for(size_t col = 1 ; col < pos_col.size(); ++col){
    std::cout << col << " " << pos_col[col]  << " " << pos_col[col-1] << std::endl;
    if(pos_col[col] - pos_col[col-1] == 2){
      // append character in between to quarks
      quarks.push_back(search.substr(pos_col[col-1]+1, 1));
    }
    else{
      // append characters to operators
      size_t dst = pos_col[col] - pos_col[col-1] - 1;
      operat.push_back(search.substr(pos_col[col-1]+1, dst));
    }
  }
  for(auto& c : quarks) std::cout << c << std::endl;
  for(auto& o : operat) std::cout << o << std::endl; 

  // now split operator strings and sort into tag
  for(size_t c = 0; c < operat.size(); ++c){

    std::vector<size_t> pos_dot;
    // get all positions of delimiters
    for(size_t pos = 0; pos < operat[c].length(); ++pos){
      char c = search[pos];
      if(c == '.') pos_dot.push_back(pos);
    }
    // sort quarks and operators
    //  loop has to end earlier, want to compare distances
//    for(size_t col = 1 ; col < pos_dot.size(); ++col){
//      std::cout << col << " " << pos_dot[col]  << " " << pos_col[col-1] << std::endl;
//      
//      if(first of splitted string =='g'){
//        append rest until dot to gamma
//      }
//
//      else if(first 1st of splitted == 'd'){
//          if(second of splitted == '0') append 000 as displacement
//          else{
//          append next three integers to displacement
//          }
//          
//      }
//
//      else if(first of splitted == p){
//          if(second of splitted == '0') append 000 as momentum
//          else{
//          append next three integers to momentum
//          }
//
//      }
//    }

  }
  // set an n array from input 
  // set_c_arr(sign., )

}





















