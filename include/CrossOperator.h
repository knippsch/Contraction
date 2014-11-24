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
                 const size_t nb, const int t_source, const int t_sink);

  void swap(const size_t nb1, const size_t nb2);

  inline const Eigen::MatrixXcd& operator()(const size_t nb, const size_t pd, 
                                            const size_t pu, const size_t dir1, 
                                            const size_t dir2, 
                                            const size_t rnd1, 
                                            const size_t rnd2, 
                                            const size_t rnd3) const {
    return X[nb][pd][pu][dir1][dir2][rnd1][rnd2][rnd3];
  }

private:
  std::vector<array_Xcd_d7_eigen> X;

};

}// end of namespace

#endif // CROSSOPERATOR_H_
