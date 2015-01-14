#ifndef TAG_HANDLING_H_
#define TAG_HANDLING_H_

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
// compose string to look for in the outfile from infile taking an input
// filename and storing the string for the search in tag
void compose_string(const char* filename, std::string& tag);

#endif //TAG_HANDLING_H_
