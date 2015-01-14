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
