#include "EigenVector.h"

void LapH::EigenVector::read_eigen_vector(const std::string& filename, const size_t t,
                                          const size_t verbose){

  const size_t dim_row = V[t].rows();
  const size_t number_of_eigen_vec = V[t].cols();

  //buffer for read in
  std::vector<std::complex<double> > eigen_vec(dim_row);

  if(verbose) std::cout << "\tReading eigenvectors from files:" << filename 
                        << std::endl;

  //setting up file
  std::ifstream infile(filename, std::ifstream::binary); 
  if (infile) {
    for (size_t nev = 0; nev < number_of_eigen_vec; ++nev) {
      infile.read( (char*) &(eigen_vec[0]), 2*dim_row*sizeof(double));
      for(size_t nrow = 0; nrow < dim_row; ++nrow){
        (V[t])(nrow, nev) = eigen_vec[nrow];
      }
    }
  }
  else {
    std::cout << "eigenvector file does not exist!!!\n" << std::endl;
    exit(0);
  }
  infile.close();

  // small test of trace and sum over the eigen vector matrix!
  if(verbose){
    std::cout << "trace of V^d*V on t = " << t << ":\t"
        << (V[t].adjoint() * V[t]).trace() << std::endl;
    std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"
        << (V[t].adjoint() * V[t]).sum() << std::endl;
  }   
}

void LapH::EigenVector::read_eigen_vector(const std::string& filename, 
                                          const size_t verbose){
  std::cout << "Reading eigenvectors!" << std::endl;
  for(int t = 0; t < V.size(); t++){
    char buff[10];
    sprintf(buff, "%03d", t);
    read_eigen_vector(filename + buff, t, verbose);
  }

}


