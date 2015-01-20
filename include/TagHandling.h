#ifndef TAG_HANDLING_H_
#define TAG_HANDLING_H_

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "IoHelpers.h"
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

// set the tag for the second message for a 2pt function given the indexpair
// of quantum numbers in op_Corr for source and sink
// void set_tag(Tag& tag, const std::pair<size_t, size_t>& i);

// set the tag for the second message for a 4pt function given the 
// indexquadruple of quantum numbers in op_Corr for source and sink
//void set_tag(Tag& tag, const std::array<size_t, 4>& i);

// compose string to look for in the outfile from infile taking an input
// filename and storing the string for the search in tag
void compose_string(const char* filename, std::string& tag);

// Write a function that turns a string to a searchable tag perhaps better to
// have this in the IoHelpers.h (Think about merge)
//
void string_to_tag(const std::string& in, Tag& out);
void tag_to_string(const Tag& in, std::string& out);
// print tag
void print_tag(const Tag& query);

// Compare two tags of correlation functions
bool compare_tags(const Tag& tag1, const Tag& tag2);

// set tags to zero n-pt functions
void zero_vec_tag(const size_t pts, std::vector<Tag> attr);
void zero_tag(const size_t pts, Tag& def);
#endif //TAG_HANDLING_H_
