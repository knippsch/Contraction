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
#include "Perambulator.h"
#include "typedefs.h"

#include "omp.h"

namespace LapH {

class Correlators{

public:
  Correlators();
  ~Correlators() {};

  void build_everything(const size_t config_i);

private:
  BasicOperator basic;
  LapH::Perambulator peram;
  std::vector<LapH::RandomVector> rnd_vec;
  LapH::VdaggerV vdaggerv;
  array_cd_d5 C4_mes;
  array_cd_d6 C2_mes;
  array_cd_d10 Corr;

  void set_corr(const size_t config){
    read_rnd_vectors_from_file(config);
    vdaggerv.build_vdaggerv(config);
    vdaggerv.build_rvdaggervr(config, rnd_vec);
    peram.read_perambulators_from_file(config);
  }
  void read_rnd_vectors_from_file (const int config_i);
  void compute_meson_corr(const int t_source, const int t_sink);
  void compute_meson_4pt_cross(LapH::CrossOperator& X, 
                               const int t_source, const int t_sink);

};

}// end of namespace

#endif // CORRELATORS_H_
