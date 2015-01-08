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
typedef std::list<std::pair<size_t, size_t> > indexlist_2;
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
    std::array<size_t,4> gamma;
    int id_VdaggerV;
    int flag_VdaggerV;
    size_t id_rVdaggerVr;
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

// struct which contains all id's from pdg which shall be used in the
// 2pt function. 
// index - pair of id's for source and sink
// The following indices could also be calculated on the fly using the id and 
// pdg. Saving them seperately is faster and easier for the output as same
// quantum numbers are often summed over.
// p_sq  - momentum squared. Identical for source/sink due to momentum 
//         conservation
// dg_so - combined dirac/gamma structure for the source
// dg_si - combined dirac/gamma structure for the sink

  struct pdg_C2 {
    size_t id;
    indexlist_2 index;
  };

// struct which contains all id's from pdg which shall be used in the
// 4pt function. 
// index - quadruple of id's for source and sink
// The following indices could also be calculated on the fly using the id and 
// pdg
// p_sq_cm - CM momentum. Identical for source/sink due to momentum. Choosing
//         this quantity != 0 corresponds to moving frames
// p_sq_so_1 - momentum squared of first source quark
// p_sq_so_2 - momentum squared of second source quark
// p_sq_si_3 - momentum squared of first sink quark
// p_sq_si_4 - momentum squared of first sink quark
// dg_so - combined dirac/gamma structure for the source
// dg_si - combined dirac/gamma structure for the sink

  struct pdg_C4 {
    size_t id;
    indexlist_4 index;
  };

typedef std::vector<pdg> vec_pdg_Corr;  
typedef std::vector<pd> vec_pd_VdaggerV;
typedef std::vector<pd_r> vec_pd_rVdaggerVr;
typedef std::vector<pdg_C2> vec_pdg_C2; 
typedef std::vector<pdg_C4> vec_pdg_C4; 

#endif // _TYPEDEFS_H_
