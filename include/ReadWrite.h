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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "boost/multi_array.hpp" 

#include "lime.h"
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
  void build_source_matrix (const int config_i, const int p_min, 
      const int p_max);
	void read_perambulators_from_file (const int config_i);
	void read_rnd_vectors_from_file (const int config_i);
  void read_lime_gauge_field_doubleprec_timeslices(const int config_i);

	Eigen::MatrixXcd* perambulator;
	Eigen::MatrixXcd*** basicoperator;
	Eigen::VectorXcd* rnd_vec;
  int* mom_squared;
  int number_of_momenta;

protected:
	std::complex<double>** momentum;

  Eigen::Matrix3cd** eigen_timeslice;
  double* gaugefield;
  int** iup;
  int** idown;
	void read_eigenvectors_from_file (Eigen::MatrixXcd& V, const int config_i, 
      const int t);
};

#endif // _READ_WRITE_H__
