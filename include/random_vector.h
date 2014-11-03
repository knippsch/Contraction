#ifndef RANDOM_VECTOR_H_
#define RANDOM_VECTOR_H_

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ranlxs.h"

namespace LapH {

typedef std::complex<double> cmplx;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class random_vector {

private:
  // the random vector
  std::vector<cmplx> vec;

public:
  // standard ctor and dtor are enough - all else is handled by std::vector
  random_vector(const size_t length) : vec(length, cmplx(0.0,0.0)) {};
  ~random_vector() {};

  // [] operator to directly access the elements of vec
  inline cmplx operator[](size_t i) const {
    return vec[i];
  }

  // computes the random vectors for the sources
  // input: seed   -> seed for the random vector
  //        length -> length of the random vector
  void set(const int seed);
  // computes the random vectors for the sources and stores them
  // input: seed   -> seed for the random vector
  //        length -> length of the random vector
  void set(const int seed, const std::string& filename);
  // ---------------------------------------------------------------------------
  //TODO: There should be functions which support the lime format!!!!!!!!!!!!!!!
  // ---------------------------------------------------------------------------
  // writing the random vector to some file
  // input: filename -> the filename
  void write_random_vector(const std::string& filename) const;
  // reading the random vector from some file
  // input: filename -> the filename
  void read_random_vector(const std::string& filename);

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif // RANDOM_VECTOR_H_ 
