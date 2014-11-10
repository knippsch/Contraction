#ifndef EIGENVECTOR_H_
#define EIGENVECTOR_H_

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense> 

namespace LapH {

typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class EigenVector {

private:
  // the eigenvectors
  vec_Xcd_eigen V;

public:
  // standard ctor and dtor are enough - all else is handled by std::vector
  EigenVector(const size_t t, const size_t row_size, const size_t col_size) : 
               V(t, Eigen::MatrixXcd(row_size, col_size)) {};
  EigenVector() {};
  ~EigenVector() {};

  // [] operator to directly access the elements of V
  inline const Eigen::MatrixXcd& operator[](size_t t) const {
    return V[t];
  }

  // ---------------------------------------------------------------------------
  //TODO: There should be functions which support the lime format!!!!!!!!!!!!!!!
  // ---------------------------------------------------------------------------
  // writing the eigen vector to some file
  // input: filename -> the filename
  // TODO: NOT IMPLEMENTED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  void write_eigen_vector(const std::string& filename) const;
  // reading the eigen vector from some file
  // input: filename -> the filename which must contain the path and the 
  //                    filename WITHOUT the timeslice (do not forget the dot :)
  void read_eigen_vector(const std::string& filename, const size_t verbose);
  // reading the eigen vector from some file for a specific timeslice
  // input: filename -> the path with the FULL filename
  //               t -> timeslice in V where the eigenvector will be written to
  void read_eigen_vector(const std::string& filename, const size_t t, 
                         const size_t verbose);

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif // EIGENVECTOR_H_ 
