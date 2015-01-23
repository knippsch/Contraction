/***************************************************************************/
//
// Short demonstration of gaugefield class
//
// Author: Christopher Helmes
//
/***************************************************************************/
#include <array>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>
#include <vector>

#include "boost/multi_array.hpp"
#include "boost/crc.hpp"
#include "CorrelatorIo2pt.h"
//#include "BasicOperator.h"
//#include "GlobalData.h"
#include "ranlxs.h"
#include "TagHandling.h"
#include "typedefs.h"

void fill_corr_rand(std::vector<cmplx>& one_corr, const int mult){
  const int T = one_corr.size();
  float re[T];
  float im[T];
    ranlxs(re, T);
    ranlxs(im, T);

    for (int t = 0; t < T; ++t){
     one_corr[t]=( std::complex<double>(re[t], im[t]) );
    }
}

int main(int ac, char* av[]) {
 size_t func = 100;
 size_t Lt = 48;
  // Global Message for the computation containing number of random vectors,
  // parameters and global checksum
 GlobalDat run_id;
 rlxs_init(2,1337);
 // set one specific tag to look for
 Tag test;
 std::string candidate;
 compose_string("tag_pars", candidate);
 string_to_tag(candidate, test);
 print_tag(test);
 // 100 Correlators
 std::vector<vec> correlators(func);
 for (auto& el : correlators) el.resize(Lt);
 std::vector<Tag> attributes(func);
 zero_vec_tag(4, attributes);
 attributes[30] = test;
 print_tag(test);
 // Fill correlators with random numbers
 for (auto& el : correlators) fill_corr_rand(el, &el-&correlators[0]);

 //write_2pt_lime("tag_check", run_id, attributes, correlators);

 //concatenate all correlation functions in one vector
 std::vector<cmplx> collect;
 for(auto& c : correlators)
   for (auto& el : c) collect.push_back(el);
 boost::crc_32_type chk_agent;
 size_t bytes = collect.size()*2*sizeof(double);
 std::cout << bytes << std::endl;
 chk_agent.process_bytes(collect.data(), bytes); 
 std::cout << "Checksum for Correlators is:" << chk_agent.checksum() << std::endl;

 std::vector<Tag> tags_in(func);
std::cout << tags_in.size() << std::endl;
 std::vector<vec> correlators_in(func);
for (auto& corr : correlators_in) corr.resize(Lt);
// std::cout << "read_in from file: big_test " << std::endl;
read_2pt_lime("tag_check", tags_in, correlators_in);
std::cout << tags_in.size() << std::endl;
 //std::vector<cmplx> result;
  //get_2pt_lime("final_write", 100, 96, id, result );
 // for(auto& dat : result) std::cout << dat << std::endl;
//  ASCII_dump_corr("checksum_test", Lt, func, 3, 2); 
  // for (auto& el : correlators_in.at(0)) std::cout << el << std::endl;
 

  return 0;
}

