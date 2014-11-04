#include "RandomVector.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::RandomVector::set(const int seed) {

  // initialisation of the rando vector to create Z2 random vector
  size_t length = vec.size();
  rlxs_init(0, seed);
  std::vector<float> rnd(2*length);
  ranlxs(&(rnd[0]), 2*length);

  // generating a Z_2 source
  for(size_t i = 0; i < length; ++i){
    const double sqrt2 = 0.5*sqrt(2.0);
    double re, im;
    if (rnd[2*i] < 0.5)
      re = -sqrt2;
    else
      re = sqrt2;
    if (rnd[2*i + 1] < 0.5)
      im = -sqrt2;
    else
      im = sqrt2;
    vec[i] = cmplx(re, im);
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::RandomVector::set(const int seed, const std::string& filename) {

  set(seed);
  write_random_vector(filename);

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::RandomVector::write_random_vector(
                                       const std::string& filename) const {
  // writing random vector to file
  FILE* fp = NULL;
  if((fp = fopen(filename.c_str(), "wb")) == NULL){
    std::cout << "failed to open file to write random vector: " 
              << filename << "\n" << std::endl;
    exit(0);
  }   
  int check_read_in = fwrite(&(vec[0]), sizeof(cmplx), vec.size(), fp);
  if(check_read_in !=  (int) vec.size())
    std::cout << "It seems that not all data are written to: "
              << filename.c_str() << "\n" << std::endl;
  fclose(fp);

} 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::RandomVector::read_random_vector(const std::string& filename) {

  // open file for reading
  FILE *fp = NULL;
  if((fp = fopen(filename.c_str(), "rb")) == NULL){
    std::cout << "failed to open file to read random vector: " 
              << filename << "\n" << std::endl;
    exit(0);
  }   
  // reading data
  int check_read_in = fread(&(vec[0]), sizeof(std::complex<double>), 
                             vec.size(), fp);
  if(check_read_in !=  (int) vec.size())
    std::cout << "It seems that not all data are written to: "
              << filename.c_str() << "\n" << std::endl;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
 
