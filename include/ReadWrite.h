/*
 * ReadWrite.h
 *
 *  Created on: Apr 13, 2014
 *      Author: werner
 */

#ifndef _READ_WRITE_H_
#define _READ_WRITE_H_

// TODO: check if they are all necessary. Doesn't matter much though
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
#include "propagator_io.h"
#include "quark.h"

/***************************Input from files**********************************/

class ReadWrite {

public:
	//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	ReadWrite ();
	virtual ~ReadWrite (); 
  void build_source_matrix ();
	void read_eigenvectors_from_file (const int config_i);
	void read_perambulators_from_file (const int config_i);
	void read_rnd_vectors_from_file (const int config_i);

	Eigen::MatrixXcd* perambulator;
	Eigen::MatrixXcd*** basicoperator;

protected:
	Eigen::VectorXcd* rnd_vec;
	Eigen::MatrixXcd* V;
	std::complex<double>** momentum;
  Eigen::MatrixXcd s;
};

#endif // _READ_WRITE_H__
