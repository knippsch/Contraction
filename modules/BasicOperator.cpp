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

static void create_gamma (struct lookup* gamma, const int i) {
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

BasicOperator::BasicOperator () {
  try{
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

    // creating gamma matrices
    gamma = new struct lookup[16];
    for(int i = 0; i < 16; ++i){
      create_gamma(gamma, i);
    }

    // memory for the perambulator, random vector and basic operator
    contraction =         new Eigen::MatrixXcd**[number_of_rnd_vec];
    contraction_dagger =  new Eigen::MatrixXcd*[number_of_rnd_vec];
    for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
      contraction[rnd_i] =        new Eigen::MatrixXcd*[number_of_rnd_vec];
      for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j) {
      contraction[rnd_i][rnd_j] = new Eigen::MatrixXcd[4];
      for(int blocknr = 0; blocknr < 4; ++blocknr){
        contraction[rnd_i][rnd_j][blocknr] =         Eigen::MatrixXcd::Zero(
            4 * number_of_eigen_vec, quarks[0].number_of_dilution_E);
      }
      }

      contraction_dagger[rnd_i] = new Eigen::MatrixXcd[4];
      for(int blocknr = 0; blocknr < 4; ++blocknr){
        // D_u^-1 = perambulator * basicoperator. Columns are always
        // the same, only permuted and multiplied with +-i or -1 by
        // gamma matrices. contraction holds the for columns
        contraction_dagger[rnd_i][blocknr] =  Eigen::MatrixXcd::Zero(
            4 * quarks[0].number_of_dilution_E, number_of_eigen_vec);
      }
    }

    // memory for (P^(b) rho V)^dagger exp(-ipx) V P^(b) rho to build u quark in 
    // charged case
    s_charged = new Eigen::MatrixXcd**[number_of_rnd_vec];
    for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
      s_charged[rnd_i] = new Eigen::MatrixXcd*[number_of_rnd_vec];
      for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j) {
        s_charged[rnd_i][rnd_j] = new Eigen::MatrixXcd[4];
        for(int blocknr = 0; blocknr < 4; blocknr++){
          s_charged[rnd_i][rnd_j][blocknr] = 
              Eigen::MatrixXcd::Zero(quarks[0].number_of_dilution_E, 
              quarks[0].number_of_dilution_E);
        }
      }
    }

    std::cout << "\t allocated memory in BasicOperator" << std::endl;

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
    delete[] contraction_dagger;

    gamma = NULL;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::~BasicOperator\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/*charged*********************************************************************/
/******************************************************************************/
/******************************************************************************/

// initializes contractions[col] with columns of D_u^-1

void BasicOperator::init_operator (const int t_source, const int t_sink, 
    ReadWrite* rewr, const char dilution, const char charge, const int p, 
    const int displ){

  clock_t t = clock();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  int t_sink_dil;

  //TODO parallelization should be possible

  switch(dilution) {

    case 'i':
      t_sink_dil = t_sink % quarks[0].number_of_dilution_T;
      break;
    case 'b':
      t_sink_dil = t_sink / quarks[0].number_of_dilution_T;
      break;
    default:
      std::cout << "Time dilution scheme not found in BasicOperator::\
        init_operator" << std::endl;
      exit(0);
    }

  // for charged particles dilute u quark V^dagger*V in rows and cols
  // s contains diluted basicoperator

  // initialize s to 0  
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j) {
      for(int blocknr = 0; blocknr < 4; ++blocknr) {
        for(int i = 0; i < quarks[0].number_of_dilution_E; i++) {
          for(int j = 0; j < quarks[0].number_of_dilution_E; j++) {
            s_charged[rnd_i][rnd_j][blocknr](i, j) = 0;
          }
        }
      }
    }
  }

  // dilute s: for dilution from right side do a left dilution with rho* 
  // and basicoperator[-p] and dagger the whole object
  // for dilution from the left, dilute again
  Eigen::MatrixXcd s_temp = Eigen::MatrixXcd::Zero(quarks[0].number_of_dilution_E,
      number_of_eigen_vec);

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
    for(int blocknr = 0; blocknr < 4; ++blocknr) {

      // s_temp holds first left dilution without daggering
      for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
        s_temp.row(vec_i % quarks[0].number_of_dilution_E) +=
            std::conj(rewr->rnd_vec[rnd_i](blocknr + vec_i * 4 + 
            4 * number_of_eigen_vec * t_sink)) *
            rewr->basicoperator[rewr->number_of_momenta - p - 1]
            [t_sink][displ].row(vec_i);
      }
      // second left dilution
      for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j) {
        for(int vec_j = 0; vec_j < number_of_eigen_vec; vec_j++) {
          s_charged[rnd_j][rnd_i][blocknr].row(vec_j % 
              quarks[0].number_of_dilution_E) +=
              std::conj(rewr->rnd_vec[rnd_j](blocknr + vec_j * 4 +
              4 * number_of_eigen_vec * t_sink)) *
              (s_temp.adjoint()).row(vec_j);
        }
      }
      // reset s_temp to 0
      for(int i = 0; i < quarks[0].number_of_dilution_E; i++) {
        for(int j = 0; j < number_of_eigen_vec; j++) {
          s_temp(i, j) = 0;
        }
      }

    }
  }
                
// works for uncharged pion, check if this can be used for charged
#if 0
  //TODO: if p == 0 skip vec_j loop (only diagonal elements)
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
  for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j){
    for(int blocknr = 0; blocknr < 4; ++blocknr) {
//    for(int blocknr2 = 0; blocknr2 < 4; ++blocknr2) {
      int blocknr2 = blocknr;
      // blocknr is identical with dirac. basicoperator blockdiagonal in 
      // diracspace -> treat every dirac index individually
      for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
        for(int vec_j = 0; vec_j < number_of_eigen_vec; ++vec_j) {
        s_charged[rnd_i][rnd_j][blocknr](vec_i % quarks[0].number_of_dilution_E,
            vec_j % quarks[0].number_of_dilution_E)
            += std::conj(rewr->rnd_vec[rnd_i](blocknr + vec_i * 4 + 
            4 * number_of_eigen_vec * t_sink)) * 
            rewr->rnd_vec[rnd_j](blocknr2 + vec_j * 4 + 4 * number_of_eigen_vec * t_sink) *
            rewr->basicoperator[p][t_sink](vec_i, vec_j);
        }
      }
  //  }
    }
  }
  }
#endif


  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {

    for(int col = 0; col < 4; ++col) {
      for(int row = 0; row  < 4; ++row){
        // propagator D_u^-1 = perambulator(tsource, tsink) * basicoperator(tsink)
        // calculate columns of D_u^-1. gamma structure can be implented by
        // reordering columns and multiplying them with constants
        for(int rnd_j = 0; rnd_j < number_of_rnd_vec; rnd_j++) {
        contraction[rnd_i][rnd_j][col].block(row * number_of_eigen_vec, 0,
            number_of_eigen_vec, quarks[0].number_of_dilution_E) =
        rewr->perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + 
            number_of_eigen_vec * row,
            (quarks[0].number_of_dilution_E) * quarks[0].number_of_dilution_D * 
            t_sink_dil + (quarks[0].number_of_dilution_E) * col,
            number_of_eigen_vec,
            (quarks[0].number_of_dilution_E)) *
        s_charged[rnd_i][rnd_j][col];
        }

        // don't need to build contraction_dagger for uncharged pion
//        if(charge != 'c') continue;

        contraction_dagger[rnd_i][row].block(col * quarks[0].number_of_dilution_E, 0,
            quarks[0].number_of_dilution_E, number_of_eigen_vec) = 
        (rewr->perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + 
            number_of_eigen_vec * row,
            (quarks[0].number_of_dilution_E) * quarks[0].number_of_dilution_D * 
            t_sink_dil + (quarks[0].number_of_dilution_E) * col,
            number_of_eigen_vec,
            (quarks[0].number_of_dilution_E))).adjoint() *
        rewr->basicoperator[rewr->number_of_momenta - p - 1][t_source][displ];
          
        // that's the best criterium I could think up for multiplication with
        // gamma_5 from left and right side. It changes the sign of the two
        // upper right and tow lower left blocks in dirac space
        if( ((row + col) == 3) || (abs(row - col) > 1) ){
        contraction_dagger[rnd_i][row].block(col * 
            quarks[0].number_of_dilution_E, 0,
            quarks[0].number_of_dilution_E, number_of_eigen_vec) *= -1;
        }
  
      }
    }

  }

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// returns D_u^-1 Gamma

void BasicOperator::get_operator_charged (Eigen::MatrixXcd**& op_1, ReadWrite* rewr, 
    const int dirac, const int t_sink){

  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = rnd_i + 1; rnd_j < number_of_rnd_vec; ++rnd_j) {

      for(int i = 0; i < 4; i++) {
        op_1[rnd_i][rnd_j].block(0, gamma[dirac].row[i] * quarks[0].number_of_dilution_E,
            4 * number_of_eigen_vec, quarks[0].number_of_dilution_E) =
        gamma[dirac].value[i] * contraction[rnd_i][rnd_j][i];
      }
    }

//check if it is faster to initialize op_1 block- instead of column-wise
#if 0
      for(int j = 0; j < 4; j++) {
        for(int i = 0; i < 4; i++){
      (op_1[rnd_i]).block(i * number_of_eigen_vec, (gamma[dirac].row[j]) * 
          number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
          gamma[dirac].value[j] * contraction[rnd_i][j + 4 * i];
      (op_1[rnd_i]).block(i * number_of_eigen_vec, (gamma[dirac].row[j]) * 
          number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
          gamma[dirac].value[j] * contraction[rnd_i][j + 4 * i];
      (op_1[rnd_i]).block(i * number_of_eigen_vec, (gamma[dirac].row[j]) * 
          number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
          gamma[dirac].value[j] * contraction[rnd_i][j + 4 * i];
      (op_1[rnd_i]).block(i * number_of_eigen_vec, (gamma[dirac].row[j]) * 
          number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
          gamma[dirac].value[j] * contraction[rnd_i][j + 4 * i];
        }
      }
#endif
    
//    }
  }

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//returns D_d^-1 Gamma

void BasicOperator::get_operator_g5 (Eigen::MatrixXcd*& op_1, const int dirac){
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

    for(int i = 0; i < 4; i++) {
      op_1[rnd_i].block(0, gamma[dirac].row[i] * number_of_eigen_vec,
//          4 * number_of_eigen_vec, number_of_eigen_vec) =
          4 * quarks[0].number_of_dilution_E, number_of_eigen_vec) =
      gamma[dirac].value[i] * contraction_dagger[rnd_i][i];
    }

//check if it is faster to initialize op_1 block- instead of column-wise
#if 0
      for(int j = 0; j < 4; j++) {
        for(int i = 0; i < 4; i++){
          (op_1[rnd_i]).block(i * number_of_eigen_vec, gamma[dirac].row[j] * 
              number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
              gamma[dirac].value[j] * contraction_dagger[rnd_i][j + 4 * i];
          (op_1[rnd_i]).block(i * number_of_eigen_vec, gamma[dirac].row[j] * 
              number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
              gamma[dirac].value[j] * contraction_dagger[rnd_i][j + 4 * i];
          (op_1[rnd_i]).block(i * number_of_eigen_vec, gamma[dirac].row[j] * 
              number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
              gamma[dirac].value[j] * contraction_dagger[rnd_i][j + 4 * i];
          (op_1[rnd_i]).block(i * number_of_eigen_vec, gamma[dirac].row[j] * 
              number_of_eigen_vec, number_of_eigen_vec, number_of_eigen_vec) = 
              gamma[dirac].value[j] * contraction_dagger[rnd_i][j + 4 * i];
        }
      }
#endif

  }

} 
 
// TODO: think about speedup from extracting factors -1 and +-i in get_operator 
// and get_operator_g5


#if 0

/******************************************************************************/
/******************************************************************************/
/*uncharged********************************************************************/
/******************************************************************************/
/******************************************************************************/

// initializes contractions[col] with columns of D_u^-1

// returns D_u^-1 Gamma

void BasicOperator::get_operator_uncharged (Eigen::MatrixXcd*& op_1, const int dirac){

  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

    for(int i = 0; i < 4; i++) {
      int rnd_j = 0;
 
      op_1[rnd_i].block(0, gamma[dirac].row[i] * number_of_eigen_vec,
          4 * number_of_eigen_vec, number_of_eigen_vec) =
      gamma[dirac].value[i] * contraction[rnd_i][rnd_j][i];
  
    }
  }
}

#endif

