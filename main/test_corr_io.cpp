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
#include <typeinfo>
#include <vector>

//#include "CorrelatorIo.h"
#include "GlobalData.h"
#include "BasicOperator.h"
#include "typedefs.h"


int main(int ac, char* av[]) {

  struct OpInfo{
    size_t id;
    std::array<int, 3> p;
    std::array<int, 3> dis_id;
    std::array<int, 4> gamma_id;

  };
  // initialize one OpInfo struct
  OpInfo op_info;
  // one three momentum
  std::array<int,3> mom {{0,1,2}};
  // one displacement vector
  std::array<int,3> der {{3,4,5}};
  // gamma configuration 
  std::array<int,4> dirac {{5,4,4,4}};
  // initialize operator info
  op_info.p = mom;
  op_info.dis_id = der;
  op_info.gamma_id = dirac;
  // print everything
  for (auto mom_comp : op_info.p) 
    std::cout << mom_comp << std::endl;
  for (auto der : op_info.dis_id) 
    std::cout << der << std::endl;
  for (auto dir : op_info.gamma_id) 
    std::cout << dir << std::endl;
  return 0;

}

