/*
 * BasicOperator.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "GlobalData.h"
#include "Perambulator.h"
#include "propagator_io.h"
#include "quark.h"
#include "typedefs.h"
#include "VdaggerV.h"

// struct for Look-up table in create_gamma and get_operator. To read as
// "in column i the row[i]-element is non-zero and its value is value[i]"
// As Gamma matrices are 4x4 matrices, row and value are 4-vectors

struct lookup {
  int row[4];
  std::complex<double> value[4];
  };

class BasicOperator {

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  BasicOperator();
  ~BasicOperator () {};

  void init_operator_u(const size_t particle_no, const size_t t_source, 
                       const char dilution, const size_t displ);
  void init_operator_d(const size_t particle_no, const size_t t_source, 
                       const char dilution, const size_t displ);
  void swap_operators();
  void get_operator_charged(array_Xcd_d2_eigen& op_1, const size_t particle_no, 
                            const size_t t_sink, const size_t dirac, 
                            const size_t p) const;
  void get_operator_g5(vec_Xcd_eigen& op_1, const size_t particle_no, 
                       const size_t t_sink, const size_t dirac, 
                       const size_t p) const;
  void get_operator_uncharged(vec_Xcd_eigen& op_1, const int particle_no, 
                              const int dirac, const int p) const;

  void read_rnd_vectors_from_file (const int config_i);

  // TODO: should be called on all functions which need that and afterward it should be freed
  //       from outside to free all the memory
  inline void set_basic(const size_t config){
    peram.read_perambulators_from_file(config);
    read_rnd_vectors_from_file(config);
    vdaggerv.build_vdaggerv(config);
    vdaggerv.build_rvdaggervr(config, rnd_vec);
  }

protected:
  LapH::Perambulator peram;
  std::vector<LapH::RandomVector> rnd_vec;
  LapH::VdaggerV vdaggerv;
  array_Xcd_d5_eigen contraction_dagger;
  array_Xcd_d6_eigen contraction;
  std::vector<struct lookup>  gamma;

};

#endif /* BASICOPERATOR_H_ */
