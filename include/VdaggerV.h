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

namespace LapH {

typedef boost::multi_array<Eigen::MatrixXcd, 3> ArrayXcdd3Eigen;

//typedef std::array<int, 3> Array;
//typedef std::vector<3dArray> MomDisVec; 
typedef std::complex<double> cmplx;
typedef boost::multi_array<cmplx, 2> ArrayCDd2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class VdaggerV {

private:
  // the onject containing [V^dagger * momentum * displacement * V]
  ArrayXcdd3Eigen vdaggerv;
  ArrayCDd2 momentum;

  void create_momenta();

public:
  VdaggerV ();
  ~VdaggerV () {};

  void build_source_matrix (const int config_i);
  // () operator to directly access the elements of vec
  inline const Eigen::MatrixXcd& operator()(const size_t p, const size_t t, 
                                      const size_t d) const {
    return vdaggerv[p][t][d];
  }

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif // VDAGGERV_H_ 


