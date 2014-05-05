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
	BasicOperator ();
	virtual ~BasicOperator ();
  void init_operator (const int t_source, const int t_sink, ReadWrite* rewr, 
      const char dilution, const int p);
  void get_operator(Eigen::MatrixXcd*& op_1, const int dirac);
  void get_operator_g5(Eigen::MatrixXcd*& op_1, const int dirac);


protected:
//  void create_gamma(struct lookup* gamma, const int dirac);
  Eigen::MatrixXcd** contraction;
  Eigen::MatrixXcd** contraction_dagger;
  struct lookup*  gamma;
};

#endif /* BASICOPERATOR_H_ */
