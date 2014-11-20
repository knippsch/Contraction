#ifndef VDAGGERV_H_
#define VDAGGERV_H_

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "EigenVector.h"
#include "GlobalData.h"
#include "RandomVector.h"

namespace LapH {

typedef boost::multi_array<Eigen::MatrixXcd, 3> ArrayXcdd3Eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> ArrayXcdd5Eigen;

typedef std::complex<double> cmplx;
typedef boost::multi_array<cmplx, 2> ArrayCDd2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class VdaggerV {

private:
  ArrayXcdd3Eigen vdaggerv;
  ArrayXcdd5Eigen rvdaggervr;
  ArrayCDd2 momentum;
  size_t nb_mom;
  bool is_vdaggerv_set;
  void create_momenta();

public:
  VdaggerV ();
  ~VdaggerV () {};

  void build_vdaggerv(const int config_i);
  void build_rvdaggervr(const int config_i, 
                        const std::vector<LapH::RandomVector>& rnd_vec);

  // return reference on vdaggerv
  inline const Eigen::MatrixXcd& return_vdaggerv(const size_t p, 
                                                 const size_t t, 
                                                 const size_t d) const {
    return vdaggerv[p][t][d];
  }
  // return reference on rvdaggervr
  inline const Eigen::MatrixXcd& return_rvdaggervr(const size_t p, 
                                                   const size_t t, 
                                                   const size_t d, 
                                                   const size_t rnd1,
                                                   const size_t rnd2) const {
    return rvdaggervr[p][t][d][rnd1][rnd2];
  }

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif // VDAGGERV_H_ 


