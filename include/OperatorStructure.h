#ifndef OPERATORSTRUCTURE_H_
#define OPERATORSTRUCTURE_H_ 

#include <iostream>
#include <iomanip>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <typeinfo>
#include <vector>

#include "GlobalData.h"
#include "typedefs.h"

#include "omp.h"

namespace Qns {

  struct pdg{
    size_t id;
    std::array<int,3> p;
    std::array<int,3> dis;
    std::array<size_t,4> gamma;
  };

  struct pdg_C2 {
    size_t p_sq;
    size_t dg_so;
    size_t dg_si;
    std::list<std::pair<size_t, size_t> > index;
  };

  struct pdg_C4 {
    size_t p_sq_so;
    size_t p_sq_si;
    //displ-gamma structure at source and at sink are coded as the same
    size_t dg_so;
    size_t dg_si;
    // generates a warning, but tuple get stressful in 
    // ../modules/CrossOperator.cpp
    std::list<std::array<size_t, 4> > index;
  };

  void init_from_infile();

  void set_Corr();
  void set_C2();
  void set_C4();
 
  void create_momenta();

  extern std::vector<pdg> op_Corr;
  extern std::vector<pdg_C2> op_C2;
  extern std::vector<pdg_C4> op_C4;

}// end of namespace

#endif // OPERATORSTRUCTURE_H_
