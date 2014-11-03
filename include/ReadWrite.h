/*
 * ReadWrite.h
 *
 *  Created on: Apr 13, 2014
 *      Author: werner
 */

#ifndef _READ_WRITE_H_
#define _READ_WRITE_H_

// TODO: check if they are all necessary. Doesn't matter much though
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/SparseCore"
//#include "boost/multi_array.hpp" 

//#include "lime.h"
#include "config_utils.h"
#include "GlobalData.h"
#include "io_utils.h"
#include "propagator_io.h"
#include "quark.h"
#include "random_vector.h"
#include "typedefs.h"

/***************************Input from files**********************************/

class ReadWrite {

friend class BasicOperator;

public:
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ReadWrite ();
  void build_source_matrix (const int config_i);
  void read_perambulators_from_file (const int config_i);
  void read_rnd_vectors_from_file (const int config_i);
//  void read_lime_gauge_field_doubleprec_timeslices(const int config_i);


protected:
  vec_Xcd_eigen perambulator;

//  vec_Xcd_eigen rnd_vec;
  std::vector<LapH::random_vector> rnd_vec;

  array_Xcd_d3_eigen basicoperator;
  array_cd_d2 momentum;
  Eigen::Matrix3cd** eigen_timeslice;
  double* gaugefield;
  int** iup;
  int** idown;

  void read_eigenvectors_from_file (Eigen::MatrixXcd& V, const int config_i, 
      const int t);
};

#endif // _READ_WRITE_H__
