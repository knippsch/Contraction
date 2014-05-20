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
#include "io_utils.h"
#include "config_utils.h"

/***************************Input from files**********************************/

class ReadWrite {

public:
	//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	ReadWrite ();
	virtual ~ReadWrite (); 
  void build_source_matrix (const int p);
	void read_eigenvectors_from_file (const int config_i);
	void read_perambulators_from_file (const int config_i);
	void read_rnd_vectors_from_file (const int config_i);
  void read_lime_gauge_field_doubleprec_timeslices(const int config_i);

	Eigen::MatrixXcd* perambulator;
	Eigen::MatrixXcd**** basicoperator;
	Eigen::MatrixXcd**** basicoperator_d;
  int number_of_momenta;

protected:
	Eigen::VectorXcd* rnd_vec;
	Eigen::MatrixXcd* V;
	Eigen::MatrixXcd* W;
	std::complex<double>** momentum;

  double* gaugefield;
  Eigen::Matrix3cd** eigen_timeslice;
  int** iup;
  int** idown;
};

#endif // _READ_WRITE_H__
