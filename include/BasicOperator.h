/*
 * BasicOperator.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "boost/multi_array.hpp"

#include "GlobalData.h"
#include "ReadWrite.h"
#include "propagator_io.h"
#include "quark.h"

typedef boost::multi_array<Eigen::MatrixXcd, 2> dim2_eigen_array;
typedef boost::multi_array<Eigen::MatrixXcd, 4> dim4_eigen_array;
typedef std::vector<Eigen::MatrixXcd> vec_eigen;
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
	BasicOperator (ReadWrite* rewr);
//	virtual ~BasicOperator ();
  void init_operator_u (const int particle_no, const int t_source, 
      const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
      const int displ);
  void init_operator_d (const int particle_no, const int t_source, 
      const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
      const int displ);
  void get_operator_charged(dim2_eigen_array& op_1, const int particle_no, const int t_sink, 
                            ReadWrite* rewr, const int dirac, const int p) const;
  void get_operator_g5(vec_eigen& op_1, const int particle_no, 
                       const int dirac, const int p) const;
  void get_operator_uncharged(vec_eigen& op_1, const int particle_no, 
                              const int dirac, const int p) const;


protected:
//  void create_gamma(struct lookup* gamma, const int dirac);
  dim4_eigen_array contraction_dagger;
  dim4_eigen_array contraction;
  std::vector<struct lookup>  gamma;
};

#endif /* BASICOPERATOR_H_ */
