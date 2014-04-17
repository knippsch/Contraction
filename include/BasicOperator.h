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

class BasicOperator {

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	BasicOperator ();
	virtual ~BasicOperator ();
  void init_operator (const int t_source, const int t_sink, ReadWrite* rewr);
  void get_operator(Eigen::MatrixXcd*& op_1);
  void get_operator_g5(Eigen::MatrixXcd*& op_1);


protected:
  Eigen::MatrixXcd** contraction;
	Eigen::SparseMatrix<std::complex<double> >* gamma;
	std::complex<double>** momentum;
};

#endif /* BASICOPERATOR_H_ */
