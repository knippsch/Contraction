#ifndef CROSSOPERATOR_H_
#define CROSSOPERATOR_H_

#include <algorithm>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <vector>

#include "GlobalData.h"
#include "BasicOperator.h"
#include "typedefs.h"
#include "VdaggerV.h"

namespace LapH {

class CrossOperator{

public:
  CrossOperator() : X() {};
  CrossOperator(const size_t number);
  ~CrossOperator() {};

  void construct(const BasicOperator& basic, const VdaggerV& vdaggerv, 
                 const size_t nb, const int t_source, const int t_sink,
                 const size_t type);

  void compute_X(const BasicOperator& basic, const size_t id_si, 
                 const Eigen::MatrixXcd& Q2, const Eigen::MatrixXcd& VdaggerV, 
                 Eigen::MatrixXcd& X);

  void swap(const size_t nb1, const size_t nb2);

  inline const Eigen::MatrixXcd& operator()(const size_t nb, 
                                            const size_t id_so,
                                            const size_t id_si, 
                                            const size_t rnd1, 
                                            const size_t rnd2, 
                                            const size_t rnd3) const {
    return X[nb][id_si][id_so][rnd1][rnd2][rnd3];
  }

private:
  std::vector<array_Xcd_d5_eigen> X;

};

}// end of namespace

#endif // CROSSOPERATOR_H_
