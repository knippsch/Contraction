#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include "GlobalData.h"
#include "OperatorStructure.h"
#include "Perambulator.h"
#include "typedefs.h"
#include "VdaggerV.h"

// struct for Look-up table in create_gamma and get_operator. To read as
// "in column i the row[i]-element is non-zero and its value is value[i]"
// As Gamma matrices are 4x4 matrices, row and value are 4-vectors

struct lookup {
  int row[4];
  std::complex<double> value[4];
};

class BasicOperator {

public:
  BasicOperator();
  ~BasicOperator () {};

  void init_operator(const char dilution, const size_t displ, 
                     const LapH::VdaggerV& vdaggerv,
                     const LapH::Perambulator& peram);

  // returns D_u^-1 Gamma D_d^-1
  inline const Eigen::MatrixXcd& get_operator(const int t1, const int t2,
                                 const int t3, const size_t index,
                                 const size_t rnd_i, const size_t rnd_j) const {
    return Q2[t1][t2][t3][index][rnd_i][rnd_j];
  }

  
  void mult_dirac(const Eigen::MatrixXcd& matrix, Eigen::MatrixXcd& reordered, 
                  const size_t index);

  std::vector<struct lookup>  gamma;

private:
  array_Xcd_d6_eigen Q2;

};

#endif /* BASICOPERATOR_H_ */
