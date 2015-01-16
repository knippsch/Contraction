#ifndef TAG_HANDLING_H_
#define TAG_HANDLING_H_

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "IoHelpers.h"
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
// compose string to look for in the outfile from infile taking an input
// filename and storing the string for the search in tag
void compose_string(const char* filename, std::string& tag);

// Write a function that turns a string to a searchable tag perhaps better to
// have this in the IoHelpers.h (Think about merge)
void string_to_tag(const std::string& search, Tag sign);
#endif //TAG_HANDLING_H_
