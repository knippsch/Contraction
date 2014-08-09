/*
 * BasicOperator.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <iostream>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "GlobalData.h"
#include "ReadWrite.h"
#include "propagator_io.h"
#include "quark.h"

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
	virtual ~BasicOperator ();
  void init_operator_u (const int particle_no, const int t_source, 
      const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
      const int displ);
  void init_operator_d (const int particle_no, const int t_source, 
      const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
      const int displ);
  void get_operator_charged(Eigen::MatrixXcd**& op_1, const int particle_no,
      const int t_sink, ReadWrite* rewr, const int dirac, const int p);
  void get_operator_g5(Eigen::MatrixXcd*& op_1, const int particle_no, 
      const int dirac, const int p);
  void get_operator_uncharged(Eigen::MatrixXcd*& op_1, const int particle_no, 
      const int dirac, const int p);


protected:
//  void create_gamma(struct lookup* gamma, const int dirac);
  Eigen::MatrixXcd**** contraction_dagger;
  Eigen::MatrixXcd**** contraction;
  struct lookup*  gamma;
  BasicOperator ();
};

#endif /* BASICOPERATOR_H_ */
