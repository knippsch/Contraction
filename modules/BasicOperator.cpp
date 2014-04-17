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

static void create_gamma (
    Eigen::SparseMatrix<std::complex<double> >* const gamma, const int i) {

  try{
    switch(i) {
    case 0: // gamm_0 (time component)
      gamma[0].insert(0, 2) = -1.;
      gamma[0].insert(1, 3) = -1.;
      gamma[0].insert(2, 0) = -1.;
      gamma[0].insert(3, 1) = -1.;
      break;
    case 1: // gamma_1
      gamma[1].insert(0, 3) = -I;
      gamma[1].insert(1, 2) = -I;
      gamma[1].insert(2, 1) = I;
      gamma[1].insert(3, 0) = I;
      break;
    case 2: // gamma_2
      gamma[2].insert(0, 3) = -1.;
      gamma[2].insert(1, 2) = 1.;
      gamma[2].insert(2, 1) = 1.;
      gamma[2].insert(3, 0) = -1.;
      break;
    case 3: // gamma_3
      gamma[3].insert(0, 2) = -I;
      gamma[3].insert(1, 3) = I;
      gamma[3].insert(2, 0) = I;
      gamma[3].insert(3, 1) = -I;
      break;
    case 4: // unity
      gamma[4].insert(0, 0) = 1.;
      gamma[4].insert(1, 1) = 1.;
      gamma[4].insert(2, 2) = 1.;
      gamma[4].insert(3, 3) = 1.;
      break;
    case 5: // gamma_5
      gamma[5].insert(0, 0) = 1.;
      gamma[5].insert(1, 1) = 1.;
      gamma[5].insert(2, 2) = -1.;
      gamma[5].insert(3, 3) = -1.;
      break;
    case 6: // gamma_0 gamma_5
      gamma[6].insert(0, 2) = 1.;
      gamma[6].insert(1, 3) = 1.;
      gamma[6].insert(2, 0) = -1.;
      gamma[6].insert(3, 1) = -1.;
      break;
    case 7: // gamma_1 gamma_5
      gamma[7].insert(0, 3) = I;
      gamma[7].insert(1, 2) = I;
      gamma[7].insert(2, 1) = I;
      gamma[7].insert(3, 0) = I;
      break;
    case 8: // gamma_2 gamma_5
      gamma[8].insert(3, 0) = 1.;
      gamma[8].insert(2, 1) = -1.;
      gamma[8].insert(1, 2) = 1.;
      gamma[8].insert(0, 3) = -1.;
      break;
    case 9: // gamma_3 gamma_5
      gamma[9].insert(0, 2) = I;
      gamma[9].insert(1, 3) = -I;
      gamma[9].insert(2, 0) = I;
      gamma[9].insert(3, 1) = -I;
      break;
    case 10: // gamma_0 gamma_1
      gamma[10].insert(0, 1) = -I;
      gamma[10].insert(1, 0) = -I;
      gamma[10].insert(2, 3) = I;
      gamma[10].insert(3, 2) = I;
      break;
    case 11: // gamma_0 gamma_2
      gamma[11].insert(0, 1) = -1.;
      gamma[11].insert(1, 0) = 1.;
      gamma[11].insert(2, 3) = 1.;
      gamma[11].insert(3, 2) = -1.;
      break;
    case 12: // gamma_0 gamma_3
      gamma[12].insert(0, 0) = -I;
      gamma[12].insert(1, 1) = I;
      gamma[12].insert(2, 2) = I;
      gamma[12].insert(3, 3) = -I;
      break;
    case 13: // gamma_1 gamma_2
      gamma[13].insert(0, 0) = I;
      gamma[13].insert(1, 1) = -I;
      gamma[13].insert(2, 2) = I;
      gamma[13].insert(3, 3) = -I;
      break;
    case 14: // gamma_1 gamma_3
      gamma[14].insert(0, 1) = -1.;
      gamma[14].insert(1, 0) = 1.;
      gamma[14].insert(2, 3) = -1.;
      gamma[14].insert(3, 2) = 1.;
      break;
    case 15: // gamma_2 gamma_3
      gamma[15].insert(0, 1) = I;
      gamma[15].insert(1, 0) = I;
      gamma[15].insert(2, 3) = I;
      gamma[15].insert(3, 2) = I;
      break;
    case 16: // gamma_2 gamma_0 gamma_5
      gamma[16].insert(0, 1) = 1.;
      gamma[16].insert(1, 0) = -1.;
      gamma[16].insert(2, 3) = 1.;
      gamma[16].insert(3, 2) = -1.;
      break;
    case 17: // nothing yet
      gamma[17].insert(0, 0) = 1.;
      gamma[17].insert(1, 1) = 1.;
      gamma[17].insert(2, 2) = 1.;
      gamma[17].insert(3, 3) = 1.;
      break;
    case 18: // nothing yet
      gamma[18].insert(0, 0) = 1.;
      gamma[18].insert(1, 1) = 1.;
      gamma[18].insert(2, 2) = 1.;
      gamma[18].insert(3, 3) = 1.;
      break;
    case 19: // nothing yet
      gamma[19].insert(0, 0) = 1.;
      gamma[19].insert(1, 1) = 1.;
      gamma[19].insert(2, 2) = 1.;
      gamma[19].insert(3, 3) = 1.;
      break;
    default:
      printf("Dirac component %d not found\n", i);
      exit(0);
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void create_momenta (std::complex<double>** momentum) {

  try{
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();
    // helper variables for momenta
    const double px = 2. * M_PI / (double) Lx;
    const double py = 2. * M_PI / (double) Ly;
    const double pz = 2. * M_PI / (double) Lz;
    int p = 0;
    // running over all momentum components
    for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
      for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
        for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
          // running over all lattice points
          for(int x = 0; x < Lx; ++x){
            const int xH = x * Ly * Lz; // helper variable
            const int ipxH = ipx * px * x; // helper variable
            for(int y = 0; y < Ly; ++y){
              const int xHyH = xH + y * Lz; // helper variable
              const int ipxHipzH = ipxH + ipy * py * y; // helper variable
              for(int z = 0; z < Lz; ++z){
                //ipz=1;
                momentum[p][xHyH + z] = exp(-I * (ipxHipzH + ipz * pz * z));
                //std::cout << "mom = " << momentum[p][xHyH + z] << std::endl;
              }
            }
          }
          ++p;
        }
      }
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_momenta\n";
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

BasicOperator::BasicOperator () {
  try{
    // TODO: Look at Lx, Ly, Lz, dim_row, number_of_max_mom, verbose
    // necassary for momenta?
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int dim_row = global_data->get_dim_row();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int number_of_max_mom = global_data->get_number_of_max_mom();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = quarks[0].number_of_dilution_T
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int verbose = global_data->get_verbose();
    // creating gamma matrices
//    gamma = new Eigen::SparseMatrix<std::complex<double> >[20];
//    for(int i = 0; i < 20; ++i){
//      (gamma[i]).resize(4, 4); // TODO: Might not be necessary. CHECK!
//      create_gamma(gamma, i);
//    }
//    if(verbose) {
//      std::cout << "gamma 0:\n" << gamma[0] << std::endl;
//      std::cout << "gamma 1:\n" << gamma[1] << std::endl;
//      std::cout << "gamma 2:\n" << gamma[2] << std::endl;
//      std::cout << "gamma 3:\n" << gamma[3] << std::endl;
//    }
    // Initializing memory for eigen vectors
    // momentum creation
    momentum = new std::complex<double>*[number_of_max_mom];
    for(int p = 0; p < number_of_max_mom; ++p)
      momentum[p] = new std::complex<double>[Lx * Ly * Lz];
    create_momenta(momentum);
    // memory for the perambulator, random vector and basic operator
    contraction = new Eigen::MatrixXcd*[number_of_rnd_vec];
    for(int i = 0; i < number_of_rnd_vec; ++i){
      contraction[i] = new Eigen::MatrixXcd[4];
      for(int col = 0; col < 4; ++col){
        // D_u^-1 = perambulator * basicoperator. Columns are always
        // the same, only permuted and multiplied with +-i or -1 by
        // gamma matrices. contraction holds the for columns
        contraction[i][col] = Eigen::MatrixXcd::Zero(
            4 * number_of_eigen_vec, number_of_eigen_vec);
      }
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::BasicOperator\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
// destructor *****************************************************************/
/******************************************************************************/
/******************************************************************************/

BasicOperator::~BasicOperator () {
  try{

    delete[] gamma;
    delete[] contraction;
    delete[] momentum;

    gamma = NULL;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::~BasicOperator\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// initializes contractions[col] with columns of D_u^-1

void BasicOperator::init_operator (const int t_source, const int t_sink, ReadWrite* rewr){

  clock_t t = clock();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  //TODO parallelization should be possible

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i)
    for(int col = 0; col < 4; ++col)
      for(int row = 0; row  < 4; ++row){
        // propagator D_u^-1 = perambulator(tsource, tsink) * basicoperator(tsink)
        // calculate columns of D_u^-1. gamma structure can be implented by
        // reordering columns and multiplying them with constants
        contraction[rnd_i][col].block(row * number_of_eigen_vec, 0, 
            number_of_eigen_vec, number_of_eigen_vec) =
        rewr->perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + 
              number_of_eigen_vec * row,
            quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D * 
            (t_sink % quarks[0].number_of_dilution_T) + 
              quarks[0].number_of_dilution_E * col,
            number_of_eigen_vec,
            quarks[0].number_of_dilution_E) *
        rewr->basicoperator[rnd_i][t_sink][col];
      }

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// returns D_u^-1 Gamma

void BasicOperator::get_operator (Eigen::MatrixXcd*& op_1){
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

  //TODO dirac has to be parameter in function header -> choose dirac structure
  //TODO implement other dirac structures

  int dirac = 5;

  switch(dirac) {
    
    case 5:
      // case gamma_5: diag(1,1,-1,-1)
      // no reordering of columns, but two factors -1
      (op_1[rnd_i]).block(0, 0,                       4 * number_of_eigen_vec, 
          number_of_eigen_vec) = contraction[rnd_i][0];
      (op_1[rnd_i]).block(0, 1 * number_of_eigen_vec, 4 * number_of_eigen_vec, 
          number_of_eigen_vec) = contraction[rnd_i][1];
      (op_1[rnd_i]).block(0, 2 * number_of_eigen_vec, 4 * number_of_eigen_vec, 
          number_of_eigen_vec) = -1 * contraction[rnd_i][2];
      (op_1[rnd_i]).block(0, 3 * number_of_eigen_vec, 4 * number_of_eigen_vec, 
          number_of_eigen_vec) = -1 * contraction[rnd_i][3];
      break;

    default:
      printf("Gamma structure %d not found\n", dirac);
      exit(0);
    } 
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//returns D_d^-1 Gamma

void BasicOperator::get_operator_g5 (Eigen::MatrixXcd*& op_1){
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

    int dirac = 5;

    switch(dirac) {
      
      case 5:
      // like in get_operator(), but with additional gamma_5 trick
      // D_d^-1 = gamma_5 D_u^-1^dagger gamma_5
      // .adjoint and filing rows rather than columns account for dagger, 
      // the changed minussigns give additional gamma_5

#if 0
      (op_1[rnd_i]).block(0                      , 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][0]).transpose();
      (op_1[rnd_i]).block(1 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][1]).transpose();
      (op_1[rnd_i]).block(2 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][2]).transpose();
      (op_1[rnd_i]).block(3 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][3]).transpose();
#endif

      (op_1[rnd_i]).block(0                      , 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][0]).adjoint();
      (op_1[rnd_i]).block(1 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][1]).adjoint();
      (op_1[rnd_i]).block(2 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][2]).adjoint();
      (op_1[rnd_i]).block(3 * number_of_eigen_vec, 0, number_of_eigen_vec, 
          4 * number_of_eigen_vec) = (contraction[rnd_i][3]).adjoint();

      (op_1[rnd_i]).block(2 * number_of_eigen_vec, 0, 2 * number_of_eigen_vec, 
          4 * number_of_eigen_vec) *= -1;
      
      break;

    default:
      printf("Gamma structure %d not found\n", dirac);
      exit(0);
    } 

  }

} 
 
// TODO: think about speedup from extracting factors -1 and +-i in get_operator 
// and get_operator_g5
