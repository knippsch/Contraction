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
#include <list>

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

// index typedefs
typedef std::list<size_t> indexlist_1;
typedef std::list<std::pair<size_t, size_t> > indexlist_2;
typedef std::list<std::array<size_t, 3> > indexlist_3;
typedef std::list<std::array<size_t, 4> > indexlist_4;

// Operator typedefs

// struct which contains all desired combinations of quantum numbers. pdg means
// momentum/displacement/gamma structure
// id    - allows to restore the physical content from the internally used id
// p3    - momentum three vector
// dis3  - displacement three vector
// gamma - gamma structure. The indices for the different gamma combinations
//         are defined in modules/BasicOperator.cpp and the used combination
//         is the product of gammas corresponding to the 4 indices defined
//         by the 4-array
// id_VdaggerV - this only runs over p/d as VdaggerV is independent of the 
//         gamma structure. Additionally, the negative momenta can be obtained 
//         from their counterpart by adjoining. To adjoin the correct momentum, 
//         the index to be used is saved seperately. 
// flag_VdaggerV - when building Q2, VdaggerV can be reused for every gamma 
//         combination as long as p/d stay the same. This is a flag which
//         is 1 for the first combination and 0 else and will be checked
//         in BasicOperator::init_operator
// id_rVdaggerVr - like id_VdaggerV, but all momenta are saved

  struct pdg{
    size_t id;
    std::array<int,3> p3;
    std::array<int,3> dis3;
    std::vector<int> gamma;
    int id_vdv;
    bool first_vdv;
    bool negative_momentum;
    size_t id_rvdvr;
  };

// struct which contains all id's from pdg for which VdaggerV is not 
// redundant. Only half the momenta and no gamma structure is contained
// id    - number to reference VdaggerV
// index - id of pdg
  struct pd{
    size_t id;
    size_t index;
  };

// struct which contains all id's from pdg for which rVdaggerVr is not 
// redundant. No gamma structure is contained. Half the momenta can be 
// calculated by adjoining their negative counterpart, which will be done if
// the adjoint-flag is set to 1
// id    - number to reference VdaggerV
// index - id of pdg and id for the adjoint pdg
  struct pd_r{
    size_t id;
    size_t id_adjoint;
    size_t index;
    bool adjoint;
  };

  struct index_2pt {
    size_t id;
    size_t index_Q2;
    size_t index_Corr;
  };

  struct index_4pt {
    size_t id;
    size_t index_Q2[2];
    size_t index_Corr[2];
  };

  struct index_IO_1 {
    size_t id;
    indexlist_1 index_pt;
  };

  struct index_IO_2 {
    size_t id;
    indexlist_2 index_pt;
  };

typedef std::vector<pdg> vec_pdg_Corr;  
typedef std::vector<pd> vec_pd_VdaggerV;
typedef std::vector<pd_r> vec_pd_rVdaggerVr;
typedef std::vector<index_2pt> vec_index_2pt;
typedef std::vector<index_4pt> vec_index_4pt;
typedef std::vector<index_IO_1> vec_index_IO_1;
typedef std::vector<index_IO_2> vec_index_IO_2;

#endif // _TYPEDEFS_H_
