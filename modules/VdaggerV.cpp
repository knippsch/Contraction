#include "VDaggerV.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
LapH::VdaggerV::VdaggerV() : vdaggerv(), momentum() {

  const int nb_ev = global_data->get_number_of_eigen_vec();
  const int nb_mom = global_data->get_number_of_momenta();
  const int Lt = global_data->get_Lt();
  const int Vs = global_data->get_Lx() * global_data->get_Ly() * 
                 global_data->get_Lz();         

  vdaggerv.resize(boost::extents[nb_mom][Lt][4]);
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.size(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  momentum.resize(boost::extents[nb_mom][Vs]);
  create_momenta();

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::create_momenta () {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();

  static const std::complex<double> I(0.0, 1.0);

  //const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();

  // helper variables for momenta
  const double px = 2. * M_PI / (double) Lx;
  const double py = 2. * M_PI / (double) Ly;
  const double pz = 2. * M_PI / (double) Lz;
  int p = 0;
  int max_mom_squared = global_data->get_number_of_max_mom();

  // running over all momentum components
  for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
    for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
      for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
        if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
          continue;
        }
        //TODO: for Lx == Ly == Lz ipxH and ipxHipyH may be integers and px, 
        //py get multiplied in the exponential
        // running over all lattice points
        for(int x = 0; x < Lx; ++x){
          const int xH = x * Ly * Lz; // helper variable
          const double ipxH = ipx * px * x; // helper variable
          for(int y = 0; y < Ly; ++y){
            const int xHyH = xH + y * Lz; // helper variable
            const double ipxHipyH = ipxH + ipy * py * y; // helper variable
            for(int z = 0; z < Lz; ++z){
              momentum[p][xHyH + z] = exp(-I * (ipxHipyH + ipz * pz * z));
            }
          }
        }
        ++p;
      }
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static void read_eigenvectors_from_file (LapH::EigenVector& V,
                                         const int config_i, const int t) {
    char name[200];
    std::string filename = global_data->get_path_eigenvectors() + "/"
        + global_data->get_name_eigenvectors();
    sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t);

    V.read_eigen_vector(name, 0, 0);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_source_matrix (const int config_i) {

  clock_t t2 = clock();
  std::cout << "\tbuild source matrix:";

  const int Lt = global_data->get_Lt();
  const int dim_row = global_data->get_dim_row();
  const int displ_min = global_data->get_displ_min();
  const int displ_max = global_data->get_displ_max();
  const int nb_ev = global_data->get_number_of_eigen_vec();
  const int nb_mom = global_data->get_number_of_momenta();

  Eigen::MatrixXcd W_t = Eigen::MatrixXcd::Zero(dim_row, nb_ev);
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);

  LapH::EigenVector V_t(1, dim_row, nb_ev);

  for(int t = 0; t < Lt; ++t){

    read_eigenvectors_from_file(V_t, config_i, t);

    for(int dir = displ_min; dir < displ_max + 1; dir++) {
      for(int p = nb_mom/2; p < nb_mom; p++){
        // TODO: checking the order of loops - enhancement might be possible
        // e.g. by reordering vdaggerv with t faster than rnd_i

        // for p = 0, s is the unit matrix
        // TODO: initialize somewhere in the constructor
        if(p == (nb_mom/2)){ // zero momentum
          (vdaggerv[nb_mom/2][t][0]) = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
        }
        else { // not zero momentum
          // momentum vector contains exp(-i p x)
          // Divisor 3 for colour index. All three colours on same lattice 
          // site get the same momentum
          for(int x = 0; x < dim_row; ++x) {
            mom(x) = momentum[p][x/3];
          }
          vdaggerv[p][t][0] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];
          // TODO: this is not needed since adjoining takes no time
          //       would save a lot of memeory
          vdaggerv[nb_mom - p - 1][t][0] = (vdaggerv[p][t][0]).adjoint();
           
         } // end if momentum
      } // end for momentum
    } // end for displacement

  } // loop over time ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
  return;
}





