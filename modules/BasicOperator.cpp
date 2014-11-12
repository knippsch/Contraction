/*
 * BasicOperator.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#include "BasicOperator.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);

// Look-up table for gamma matrices. For every Gamma structure (currently 0-15)
// the four non-zero values are specified.

static void create_gamma (std::vector<struct lookup>& gamma, const int i) {
  try {
    switch(i) {
    case 0: // gamma_0
      gamma[0].row[0] = 2;
      gamma[0].value[0] = -1;
      gamma[0].row[1] = 3;
      gamma[0].value[1] = -1;
      gamma[0].row[2] = 0;
      gamma[0].value[2] = -1;
      gamma[0].row[3] = 1;
      gamma[0].value[3] = -1;
      break;

    case 1: // gamma_1
      gamma[1].row[0] = 3;
      gamma[1].value[0] = I;
      gamma[1].row[1] = 2;
      gamma[1].value[1] = I;
      gamma[1].row[2] = 1;
      gamma[1].value[2] = -I;
      gamma[1].row[3] = 0;
      gamma[1].value[3] = -I;
      break;

    case 2: // gamma_2
      gamma[2].row[0] = 3;
      gamma[2].value[0] = -1;
      gamma[2].row[1] = 2;
      gamma[2].value[1] = 1;
      gamma[2].row[2] = 1;
      gamma[2].value[2] = 1;
      gamma[2].row[3] = 0;
      gamma[2].value[3] = -1;
      break;

    case 3: // gamma_3
      gamma[3].row[0] = 2;
      gamma[3].value[0] = I;
      gamma[3].row[1] = 3;
      gamma[3].value[1] = -I;
      gamma[3].row[2] = 0;
      gamma[3].value[2] = -I;
      gamma[3].row[3] = 1;
      gamma[3].value[3] = I;
      break;

    case 4: // unity
      gamma[4].row[0] = 0;
      gamma[4].value[0] = 1;
      gamma[4].row[1] = 1;
      gamma[4].value[1] = 1;
      gamma[4].row[2] = 2;
      gamma[4].value[2] = 1;
      gamma[4].row[3] = 3;
      gamma[4].value[3] = 1;
      break;

    case 5: // gamma_5
      gamma[5].row[0] = 0;
      gamma[5].value[0] = 1;
      gamma[5].row[1] = 1;
      gamma[5].value[1] = 1;
      gamma[5].row[2] = 2;
      gamma[5].value[2] = -1;
      gamma[5].row[3] = 3;
      gamma[5].value[3] = -1;
      break;

    case 6: // gamma_0 * gamma_5
      gamma[6].row[0] = 2;
      gamma[6].value[0] = -1;
      gamma[6].row[1] = 3;
      gamma[6].value[1] = -1;
      gamma[6].row[2] = 0;
      gamma[6].value[2] = 1;
      gamma[6].row[3] = 1;
      gamma[6].value[3] = 1;
      break;

    case 7: // gamma_1 * gamma_5
      gamma[7].row[0] = 3;
      gamma[7].value[0] = I;
      gamma[7].row[1] = 2;
      gamma[7].value[1] = I;
      gamma[7].row[2] = 1;
      gamma[7].value[2] = I;
      gamma[7].row[3] = 0;
      gamma[7].value[3] = I;
      break;

    case 8: // gamma_2 * gamma_5
      gamma[8].row[0] = 3;
      gamma[8].value[0] = -1;
      gamma[8].row[1] = 2;
      gamma[8].value[1] = 1;
      gamma[8].row[2] = 1;
      gamma[8].value[2] = -1;
      gamma[8].row[3] = 0;
      gamma[8].value[3] = 1;
      break;

    case 9: // gamma_3 * gamma_5
      gamma[9].row[0] = 2;
      gamma[9].value[0] = I;
      gamma[9].row[1] = 3;
      gamma[9].value[1] = -I;
      gamma[9].row[2] = 0;
      gamma[9].value[2] = I;
      gamma[9].row[3] = 1;
      gamma[9].value[3] = -I;
      break;

    case 10: // gamma_0 * gamma_1
      gamma[10].row[0] = 1;
      gamma[10].value[0] = -I;
      gamma[10].row[1] = 0;
      gamma[10].value[1] = -I;
      gamma[10].row[2] = 3;
      gamma[10].value[2] = I;
      gamma[10].row[3] = 2;
      gamma[10].value[3] = I;
      break;

    case 11: // gamma_0 * gamma_2
      gamma[11].row[0] = 1;
      gamma[11].value[0] = 1;
      gamma[11].row[1] = 0;
      gamma[11].value[1] = -1;
      gamma[11].row[2] = 3;
      gamma[11].value[2] = -1;
      gamma[11].row[3] = 2;
      gamma[11].value[3] = 1;
      break;

    case 12: // gamma_0 * gamma_3
      gamma[12].row[0] = 0;
      gamma[12].value[0] = -I;
      gamma[12].row[1] = 1;
      gamma[12].value[1] = I;
      gamma[12].row[2] = 2;
      gamma[12].value[2] = I;
      gamma[12].row[3] = 3;
      gamma[12].value[3] = -I;
      break;

    case 13: // gamma_1 * gamma_2
      gamma[13].row[0] = 0;
      gamma[13].value[0] = I;
      gamma[13].row[1] = 1;
      gamma[13].value[1] = -I;
      gamma[13].row[2] = 2;
      gamma[13].value[2] = I;
      gamma[13].row[3] = 3;
      gamma[13].value[3] = -I;
      break;

    case 14: // gamma_1 * gamma_3
      gamma[14].row[0] = 1;
      gamma[14].value[0] = 1;
      gamma[14].row[1] = 0;
      gamma[14].value[1] = -1;
      gamma[14].row[2] = 3;
      gamma[14].value[2] = 1;
      gamma[14].row[3] = 2;
      gamma[14].value[3] = -1;
      break;

    case 15: // gamma_2 * gamma_3
      gamma[15].row[0] = 1;
      gamma[15].value[0] = I;
      gamma[15].row[1] = 0;
      gamma[15].value[1] = I;
      gamma[15].row[2] = 3;
      gamma[15].value[2] = I;
      gamma[15].row[3] = 2;
      gamma[15].value[3] = I;
      break;

    case 16: // gamma_2 * gamma_0 * gamma_5
      gamma[16].row[0] = 1;
      gamma[16].value[0] = -1;
      gamma[16].row[1] = 0;
      gamma[16].value[1] = 1;
      gamma[16].row[2] = 3;
      gamma[16].value[2] = -1;
      gamma[16].row[3] = 2;
      gamma[16].value[3] = 1;
      break;
    default:
      printf("Dirac component %d not found in BasicOperator::create_gamma\n", i);
      exit(0);
    }
  return;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
} // internal namespace ends here

/******************************************************************************/
/******************************************************************************/
// constructor ****************************************************************/
/******************************************************************************/
/******************************************************************************/

BasicOperator::BasicOperator() : peram(),
                                 rnd_vec(),
                                 vdaggerv(),
                                 contraction_dagger(), 
                                 contraction(), 
                                 gamma() {

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;

  // creating gamma matrices
  gamma.resize(16);
  for(int i = 0; i < 16; ++i)
    create_gamma(gamma, i);

  // D_u^-1 = perambulator * basicoperator. Columns are always
  // the same, only permuted and multiplied with +-i or -1 by
  // gamma matrices. contraction holds the for columns, contraction_dagger
  // holds the columns after gamma_5 trick
  // TODO: Last index will be dirac index and must be changed in the future
  contraction.resize(boost::extents[2][Lt][nb_mom][nb_rnd][nb_rnd][1]);
  std::fill(contraction.data(), contraction.data() + 
                                contraction.num_elements(), 
                                Eigen::MatrixXcd::Zero(4 * nb_ev, 4 * dilE));
  contraction_dagger.resize(boost::extents[2][Lt][nb_mom][nb_rnd][1]);
  std::fill(contraction_dagger.data(), contraction_dagger.data() + 
                                       contraction_dagger.num_elements(), 
                                       Eigen::MatrixXcd::Zero(4 * dilE, 4 * nb_ev));

  // TODO: the resize is unnasassarry if it is done in initialiser list 
  // with the ctor
  rnd_vec.resize(nb_rnd, LapH::RandomVector(Lt*nb_ev*4));

  std::cout << "\tallocated memory in BasicOperator" << std::endl;

}

/******************************************************************************/
/******************************************************************************/
// destructor *****************************************************************/
/******************************************************************************/
/******************************************************************************/

// initializes contractions[col] with columns of D_u^-1

void BasicOperator::init_operator_u (const size_t particle_no, 
                                     const size_t t_source, 
                                     const char dilution, const size_t displ){

  const size_t Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t nb_mom = global_data->get_number_of_momenta();

  #pragma omp parallel 
  {
  Eigen::MatrixXcd M(4 * nb_ev, 4 * dilE);
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){
  // TODO: just a workaround
  size_t t_sink_dil;
  switch(dilution) {

    case 'i':
      t_sink_dil = t % quarks[0].number_of_dilution_T;
      break;
    case 'b':
      t_sink_dil = t / quarks[0].number_of_dilution_T;
      break;
    default:
      std::cout << "Time dilution scheme not found in BasicOperator::\
        init_operator" << std::endl;
      exit(0);
  }
  // TODO: think about dirac structure and the last index - we do not want
  //       to calculate these objects too often
  // TODO: Dirac loop must vanish 
  for(size_t p = 0; p < nb_mom; p++){
    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
    for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) { 
    if(rnd_i != rnd_j){
      for(size_t col = 0; col < 4; ++col) {
      for(size_t row = 0; row  < 4; ++row){ 
        // calculate columns of D_u^-1. gamma structure can be implented by
        // reordering columns and multiplying them with constants 
        if(p <= nb_mom/2){
          M.block(row * nb_ev, col * dilE, nb_ev, dilE) =
            peram[rnd_i].block((4 * t_source + row) * nb_ev,
                               dilE * (quarks[0].number_of_dilution_D * 
                               t_sink_dil + col), nb_ev, dilE) * 
             (vdaggerv.return_rvdaggervr(p, t, 0, rnd_i, rnd_j))
                              .block(0, col*dilE, dilE, dilE);
        }
        else {
          M.block(row * nb_ev, col * dilE, nb_ev, dilE) =
            peram[rnd_i].block((4 * t_source + row) * nb_ev,
                               dilE * (quarks[0].number_of_dilution_D * 
                               t_sink_dil + col), nb_ev, dilE) * 
             (vdaggerv.return_rvdaggervr(nb_mom-p-1, t, 0, rnd_j, rnd_i)
                              .block(0, col*dilE, dilE, dilE)).adjoint();
        }
      }}// loops over col and row
      for(size_t dirac = 5; dirac < 6; dirac++){    
      for(size_t col = 0; col < 4; col++) {
        contraction[particle_no][t][p][rnd_i][rnd_j][0]
                          .block(0, col*dilE, 4*nb_ev, dilE) =
          gamma[dirac].value[col] * 
          M.block(0, gamma[dirac].row[col]*dilE, 4*nb_ev, dilE);
      }}// loop dirac ends here
    }}}// loops over rnd_i and rnd_j
  }}// loops over momenta and time end here
  }// pragma omp ends
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// initializes contractions_dagger[col] with columns of D_d^-1
void BasicOperator::init_operator_d (const size_t particle_no, 
                                     const size_t t_source, 
                                     const char dilution, const size_t displ){

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t nb_mom = global_data->get_number_of_momenta();

  #pragma omp parallel 
  {
  Eigen::MatrixXcd M(4 * dilE, 4 * nb_ev);
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){
  // TODO: just a workaround
  int t_sink_dil;
  switch(dilution) {

    case 'i':
      t_sink_dil = t % quarks[0].number_of_dilution_T;
      break;
    case 'b':
      t_sink_dil = t / quarks[0].number_of_dilution_T;
      break;
    default:
      std::cout << "Time dilution scheme not found in BasicOperator::\
        init_operator" << std::endl;
      exit(0);
  }
  // TODO: think about dirac structure and the last index - we do not want
  //       to calculate these objects too often
  // TODO: Dirac loop must vanish 
  for(size_t p = 0; p < nb_mom; p++){
    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(size_t col = 0; col < 4; ++col) {
      for(size_t row = 0; row  < 4; ++row){
        // propagator D_d^-1 = perambulator(t_source, t_sink)^dagger * 
        // basicoperator(tsource) 
        // = gamma_5 D_u^-1 gamma_5 according to gamma_5 trick
        // only necassary to build this for charged particles.
        // TODO: implement a flag to omit this calculation
        // TODO: implement more versatile momentum structure
        if(p <= nb_mom/2){
          M.block(col * dilE, nb_ev * row, dilE, nb_ev) = 
            (peram[rnd_i].block(4 * nb_ev * t_source + nb_ev * row, 
                      dilE * quarks[0].number_of_dilution_D * t_sink_dil + 
                      dilE * col, nb_ev, dilE)).adjoint() *
            vdaggerv.return_vdaggerv(p, t_source, displ);
        }
        else {
          M.block(col * dilE, row * nb_ev, dilE, nb_ev) = 
            (peram[rnd_i].block(4 * nb_ev * t_source + nb_ev * row,
                      dilE * quarks[0].number_of_dilution_D * t_sink_dil + 
                      dilE * col, nb_ev, dilE)).adjoint() *
            (vdaggerv.return_vdaggerv(nb_mom-p-1, t_source, displ)).adjoint();
        }  
        // gamma_5 trick. It changes the sign of the two upper right and two
        // lower left blocks in dirac space
        if( ((row + col) == 3) || (abs(row - col) > 1) )
          M.block(col * dilE, row * nb_ev, dilE, nb_ev) *= -1.;
      }}// loops over row and col end here
      for(size_t dirac = 5; dirac < 6; dirac++){    
      for(size_t col = 0; col < 4; col++) {
        contraction_dagger[particle_no][t][p][rnd_i][0]
                          .block(0, col*nb_ev, 4*dilE, nb_ev) =
          gamma[dirac].value[col] * 
          M.block(0, gamma[dirac].row[col]*nb_ev, 4*dilE, nb_ev);
      }}// loop dirac ends here
    }// loop over rnd_i ends here
  }}// loops over momenta and time end here
  }// pragma omp ends
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::swap_operators(){
  
  const size_t Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t nb_mom = global_data->get_number_of_momenta();
  // TODO: change the col entry to dirac
  #pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){
  for(size_t p = 0; p < nb_mom; p++){
    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(size_t col = 0; col < 1; ++col) {
        contraction_dagger[0][t][p][rnd_i][col].swap(
        contraction_dagger[1][t][p][rnd_i][col]);
      }
      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) { 
        for(size_t col = 0; col < 1; ++col) {
          contraction[0][t][p][rnd_i][rnd_j][col].swap(
          contraction[1][t][p][rnd_i][rnd_j][col]);
        }
      }
    }
  }}
}
/******************************************************************************/
/******************************************************************************/
/*uncharged********************************************************************/
/******************************************************************************/
/******************************************************************************/

// returns D_u^-1 Gamma
void BasicOperator::get_operator_uncharged (vec_Xcd_eigen& op_1, 
    const int particle_no, const int dirac, const int p) const{

//  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
//  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
//
//  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
//  for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j){
//
//    for(int i = 0; i < 4; i++) {
// 
//      op_1[rnd_i].block(0, gamma[dirac].row[i] * number_of_eigen_vec,
//          4 * number_of_eigen_vec, number_of_eigen_vec) =
//      gamma[dirac].value[i] * contraction[particle_no][p][rnd_i][rnd_j][i];
//  
//    }}
//  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::read_rnd_vectors_from_file (const int config_i) {

  clock_t t = clock();
  char infile[400];
  const int Lt = global_data->get_Lt();
  const int verbose = global_data->get_verbose();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int rnd_vec_length = Lt * number_of_eigen_vec * 4;

  char temp[100];

  if(verbose){
    std::cout << "\treading randomvectors from files:" << std::endl;
  } else {
    std::cout << "\treading randomvectors:";
  }

  int check_read_in = 0;
  for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
    // data path Christians perambulators
//      std::string filename = global_data->get_path_perambulators() + "/";

    // data path for qbig contractions
    sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
    std::string filename = global_data->get_path_perambulators()
      + "/" + temp;

    // data path for juqueen contractions
//      sprintf(temp, "cnfg%d/", config_i);
//      std::string filename = global_data->get_path_perambulators()
//				+ "/" + temp;

    // read random vector
    sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev%04d.%04d", 
        filename.c_str(), rnd_vec_i, number_of_eigen_vec, config_i);
//      sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev0120.%04d", 
//          filename.c_str(), rnd_vec_i, config_i);

//      sprintf(infile, "%s.%03d.u.Ti.%04d", filename.c_str(), rnd_vec_i,
//          config_i);

    // TODO:: explicit type conversion - Bad style
    rnd_vec[rnd_vec_i].read_random_vector(infile);
  }
  t = clock() - t;
  if(!verbose) std::cout << "\t\tSUCCESS - " << std::fixed 
                         << std::setprecision(1)
                         << ((float) t)/CLOCKS_PER_SEC << " seconds" 
                         << std::endl; 
}


