/*
 * typedefs.h
 *
 *  Created on: Aug 28, 2014
 *      Author: C. Jost
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <algorithm>
#include <complex>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "boost/multi_array.hpp"

// cpp standard library typedefs
typedef std::complex<double> cmplx;
typedef std::vector<cmplx> vec;
typedef boost::multi_array<cmplx, 2> array_cd_d2;
typedef boost::multi_array<cmplx, 3> array_cd_d3;
typedef boost::multi_array<cmplx, 4> array_cd_d4;
typedef boost::multi_array<cmplx, 5> array_cd_d5;
typedef boost::multi_array<cmplx, 6> array_cd_d6;
typedef boost::multi_array<cmplx, 7> array_cd_d7;
typedef boost::multi_array<cmplx, 8> array_cd_d8;
typedef boost::multi_array<cmplx, 9> array_cd_d9;
typedef boost::multi_array<cmplx, 10> array_cd_d10;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_Xcd_d2_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 3> array_Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 4> array_Xcd_d4_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> array_Xcd_d5_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 6> array_Xcd_d6_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 7> array_Xcd_d7_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 8> array_Xcd_d8_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 9> array_Xcd_d9_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 10> array_Xcd_d10_eigen;

// Operator typedefs
  struct pdg{
    size_t id;
    std::array<int,3> p;
    std::array<int,3> dis;
    std::array<size_t,4> gamma;
  };

  struct pdg_C2 {
    size_t p_sq;
    size_t dg_so;
    size_t dg_si;
    std::list<std::pair<size_t, size_t> > index;
  };

  struct pdg_C4 {
    size_t p_sq_so;
    size_t p_sq_si;
    //displ-gamma structure at source and at sink are coded as the same
    size_t dg_so;
    size_t dg_si;
    std::list<std::array<size_t, 4> > index;
  };


typedef std::vector<pdg> vec_pdg_Corr;  
typedef std::vector<pdg_C2> vec_pdg_C2; 
typedef std::vector<pdg_C4> vec_pdg_C4; 

#endif // _TYPEDEFS_H_
