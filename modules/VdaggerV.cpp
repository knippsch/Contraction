#include "VdaggerV.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
namespace {

inline void read_eigenvectors_from_file(LapH::EigenVector& V, 
                                        const int config_i, const int t) {
  char name[200];
  std::string filename = global_data->get_path_eigenvectors() + "/"
      + global_data->get_name_eigenvectors();
  sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t);

  V.read_eigen_vector(name, 0, 0);
}

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

LapH::VdaggerV::VdaggerV() : vdaggerv(), rvdaggervr(), momentum(),
                             is_vdaggerv_set(false) {

  const size_t Lt = global_data->get_Lt();
  const size_t Vs = global_data->get_Lx() * global_data->get_Ly() * 
                    global_data->get_Lz();         
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  const vec_pd_VdaggerV op_VdaggerV = global_data->get_lookup_VdaggerV();
  const vec_pd_rVdaggerVr op_rVdaggerVr = global_data->get_lookup_rVdaggerVr();

  const size_t nb_VdaggerV = op_VdaggerV.size();
  const size_t nb_rVdaggerVr = op_rVdaggerVr.size();

  // only half of the array is stored to save memory. But be careful, it 
  // must be mapped correctly from outside by addressing the momentum
  // correctly and daggering
  vdaggerv.resize(boost::extents[nb_VdaggerV][Lt]);
  rvdaggervr.resize(boost::extents[nb_rVdaggerVr][Lt][nb_rnd][nb_rnd]);

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(boost::extents[op_VdaggerV.size()][Vs]);
  create_momenta();

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::create_momenta () {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();

  const vec_pd_VdaggerV op_VdaggerV = global_data->get_lookup_VdaggerV();
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  static const std::complex<double> I(0.0, 1.0);

  // To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
  // quantum number id in op_VdaggerV will be used. The rest can be obtained
  // by adjoining
  for(const auto& op : op_VdaggerV){

    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op_Corr[op.index].p3[0] * 2. * M_PI / (double) Lx; 
    const double ipy = op_Corr[op.index].p3[1] * 2. * M_PI / (double) Ly;
    const double ipz = op_Corr[op.index].p3[2] * 2. * M_PI / (double) Lz;

    // calculate \vec{p} \cdot \vec{x} for all \vec{x} on the lattice
    for(int x = 0; x < Lx; ++x){
      const int xH = x * Ly * Lz; // helper variable
      const double ipxH = ipx * x; // helper variable
      for(int y = 0; y < Ly; ++y){
        const int xHyH = xH + y * Lz; // helper variable
        const double ipxHipyH = ipxH + ipy * y; // helper variable
        for(int z = 0; z < Lz; ++z){
          // multiply \vec{p} \cdot \vec{x} with complex unit and exponentiate
          momentum[op.id][xHyH + z] = exp(-I * (ipxHipyH + ipz * z));
    }}}//loops over spatial vectors end here
  }//loop over redundant quantum numbers ends here

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_vdaggerv (const int config_i) {

  clock_t t2 = clock();
  std::cout << "\tbuild vdaggerv:";

  const size_t Lt = global_data->get_Lt();
  const size_t dim_row = global_data->get_dim_row();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t id_unity = global_data->get_index_of_unity();

  const vec_pd_VdaggerV op_VdaggerV = global_data->get_lookup_VdaggerV();

  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    read_eigenvectors_from_file(V_t, config_i, t);

    // VdaggerV is independent of the gamma structure and momenta connected by
    // sign flip are related by adjoining VdaggerV. Thus the expensive 
    // calculation must only be performed for a subset of quantum numbers given
    // in op_VdaggerV.
    for(const auto& op : op_VdaggerV){

      // For zero momentum and displacement VdaggerV is the unit matrix, thus
      // the calculation is not performed
      if(op.index != id_unity){
        // momentum vector contains exp(-i p x). Divisor 3 for colour index. 
        // All three colours on same lattice site get the same momentum.
        for(size_t x = 0; x < dim_row; ++x) {
          mom(x) = momentum[op.id][x/3];
        }
        vdaggerv[op.id][t] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];

      }
      else{
        // zero momentum 
        (vdaggerv[op.id][t]) = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }

    }

//    }} // loop over momentum and displacement
  } // loop over time
}// pragma omp parallel ends here

  // set flag that vdaggerv is set
  is_vdaggerv_set = true;


  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_rvdaggervr(const int config_i,
                               const std::vector<LapH::RandomVector>& rnd_vec) {

  // check of vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
              << " computed\n\n" << std::endl;
    exit(0);
  }

  clock_t t2 = clock();
  std::cout << "\tbuild rvdaggervr:";

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  const vec_pd_rVdaggerVr op_rVdaggerVr = global_data->get_lookup_rVdaggerVr();
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  std::fill(rvdaggervr.data(), rvdaggervr.data() + rvdaggervr.num_elements(), 
            Eigen::MatrixXcd::Zero(dilE, 4*dilE));

  // TODO: just a workaround
  // can be changed to op by running over p = op/nb_dg, but dis currently
  // not supported.

  #pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){

  // rvdaggervr is calculated by multiplying the vdaggerv with the same quantum
  // numbers with random vectors from right and left. rvdaggervr for momenta 
  // with opposing sign are related by adjoining. Thus it suffices to calculate
  // it for half of the indices for which the flag op.adjoint < 0.
  for(const auto& op : op_rVdaggerVr){
    if(op.adjoint == false){

      size_t id_VdaggerV = op_Corr[op.index].id_vdv;

      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nb_ev, 4*dilE);
        // dilution from left
        for(size_t block= 0; block < 4; block++){
        for(size_t vec_i = 0; vec_i < nb_ev; ++vec_i) {
          size_t blk_i =  block + vec_i * 4 + 4 * nb_ev * t;
          
          M.block(0, vec_i%dilE + dilE*block, nb_ev, 1) += 
               vdaggerv[id_VdaggerV][t].col(vec_i) * 
               rnd_vec[rnd_i][blk_i];
        }}// end of dilution
        for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
        if(rnd_i != rnd_j){
          // dilution from right
          for(size_t block = 0; block < 4; block++){
          for(size_t vec_j = 0; vec_j < nb_ev; ++vec_j) {
            size_t blk_j =  block + vec_j * 4 + 4 * nb_ev * t;
            rvdaggervr[op.id][t][rnd_j][rnd_i]
                          .block(vec_j%dilE, dilE*block , 1, dilE) +=
                M.block(vec_j, dilE*block, 1, dilE) * 
                std::conj(rnd_vec[rnd_j][blk_j]);
          }}// end of dilution
        }}// rnd_j loop ends here
      }// rnd_i loop ends here
    }
  }

  // rvdaggervr for momenta with opposing sign are related by adjoining. Thus
  // for half of the indices, the calculation reduces to adjoining the
  // corresponding rvdaggervr and swapping the random vectors (as the order
  // of multiplication is reversed). The index of corresponding quantum numbers
  // is op.adjoint. It serves as flag for adjoining simultaneously, as it is
  // positive if and only if it shall be adjoined.
  // Need to loop twice as the corresponding rvdaggervr must all be calculated
  // already.
  for(const auto& op : op_rVdaggerVr){
    if(op.adjoint == true){

      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
      if(rnd_i != rnd_j){

        // rvdaggervr is a blockdiagonal 4*dilE x 4*dilE matrix. To save memory,
        // only the diagonal blocks are saved and it is written as a column 
        // vector of blocks. To reproduce the correct behavior under adjoining, 
        // the blocks have to be adjoined seperately.
        // is .adjoint().transpose() faster?
        for(size_t block = 0; block < 4; block++){
          rvdaggervr[op.id][t][rnd_j][rnd_i]
                              .block(0, block*dilE, dilE, dilE) =
            (rvdaggervr[op.id_adjoint][t][rnd_i][rnd_j]
                              .block(0, block*dilE, dilE, dilE)).adjoint();
        }
      }}}// loops over rnd vecs

    }
  }

  }// time, momemtum and displacement loop ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;

}

