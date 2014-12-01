#ifndef CORRELATORS_H_
#define CORRELATORS_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "CrossOperator.h" 
#include "GlobalData.h"
#include "BasicOperator.h"
#include "OperatorStructure.h"
#include "Perambulator.h"
#include "typedefs.h"

#include "omp.h"

namespace LapH {

class Correlators{

public:
  Correlators();
  ~Correlators() {};

  void compute_correlators(const size_t config_i);

private:
  BasicOperator basic;
  LapH::Perambulator peram;
  std::vector<LapH::RandomVector> rnd_vec;
  LapH::VdaggerV vdaggerv;
  array_cd_d5 C4_mes;
  array_cd_d4 C2_mes;
  array_cd_d6 Corr;

  void set_corr(const size_t config){
    read_rnd_vectors_from_file(config);
    vdaggerv.build_vdaggerv(config);
    vdaggerv.build_rvdaggervr(config, rnd_vec);
    peram.read_perambulators_from_file(config);
  }
  void read_rnd_vectors_from_file (const int config_i);
  void compute_meson_small_traces(const int t_source, const int t_sink);
  void compute_meson_4pt_cross_trace(LapH::CrossOperator& X, 
                                     const int t_source, const int t_sink);
  void build_and_write_2pt(const size_t config_i);
  void write_C4_3(const size_t config_i);
  void build_and_write_C4_1(const size_t config_i);
  void build_and_write_C4_2(const size_t config_i);

};

}// end of namespace

#endif // CORRELATORS_H_
