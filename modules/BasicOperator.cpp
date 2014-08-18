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

BasicOperator::BasicOperator() : contraction_dagger(), contraction(), gamma() {
  try{
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_momenta = global_data->get_number_of_momenta();

    // creating gamma matrices
    gamma.resize(16);
    for(int i = 0; i < 16; ++i)
      create_gamma(gamma, i);

    // memory for the perambulator, random vector and basic operator
    // D_u^-1 = perambulator * basicoperator. Columns are always
    // the same, only permuted and multiplied with +-i or -1 by
    // gamma matrices. contraction holds the for columns, contraction_dagger
    // holds the columns after gamma_5 trick

    const size_t nmom = number_of_momenta;
    const size_t nrnd = number_of_rnd_vec;
    contraction.resize(boost::extents[2][nmom][nrnd][4]);
    contraction_dagger.resize(boost::extents[2][nmom][nrnd][4]);

    for(int particle_no = 0; particle_no < 2; particle_no++){
      for(int p = 0; p < number_of_momenta; p++){
        for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
          for(int blocknr = 0; blocknr < 4; ++blocknr){
            // for charged pion the necessary matrix size is 
            // 4 * quarks[0].number_of_dilution_E, number_of_eigen_vec
            // but to use the same object as for pi0, the allocated size 
            // is larger. Could unite contraction_dagger and contraction
            contraction[particle_no][p][rnd_i][blocknr] = 
                Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec, 
                number_of_eigen_vec);
            contraction_dagger[particle_no][p][rnd_i][blocknr] =  
                Eigen::MatrixXcd::Zero(4 * quarks[0].number_of_dilution_E, 
                number_of_eigen_vec);
          }
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

//BasicOperator::~BasicOperator () {
//  try{
//    const std::vector<quark> quarks = global_data->get_quarks();
//    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
//    const int number_of_momenta = global_data->get_number_of_momenta();
//    delete[] gamma;
//    for(int particle_no = 0; particle_no < 2; particle_no++){
//      for(int p = 0; p < number_of_momenta; p++){
//        for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
//          delete [] contraction[particle_no][p][rnd_i];
//          delete [] contraction_dagger[particle_no][p][rnd_i];
//        }
//        delete [] contraction[particle_no][p];
//        delete [] contraction_dagger[particle_no][p];
//      }
//      delete [] contraction[particle_no];
//      delete [] contraction_dagger[particle_no];
//    }
//    delete [] contraction;
//    delete [] contraction_dagger;
//  }
//  catch(std::exception& e){
//    std::cout << e.what() << "in: BasicOperator::~BasicOperator\n";
//    exit(0);
//  }
//}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// initializes contractions[col] with columns of D_u^-1

void BasicOperator::init_operator_u (const int particle_no, const int t_source, 
    const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
    const int displ){

  clock_t t = clock();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const vec_eigen rnd_vec = rewr->get_random_vector(); 
  const dim3_eigen_array basicoperator = rewr->get_basicoperator(); 
  const vec_eigen perambulator = rewr->get_perambulator(); 

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
  // memory for (P^(b) rho V)^dagger exp(-ipx) V P^(b) rho to build u quark in 
  // s contains diluted basicoperator

  dim2_eigen_array s(boost::extents[number_of_rnd_vec][4]);
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int blocknr = 0; blocknr < 4; blocknr++){
      s[rnd_i][blocknr] = 
          Eigen::MatrixXcd::Zero(quarks[0].number_of_dilution_E, number_of_eigen_vec);
    }
  }

  // TODO: put dilution into main loop

//  // initialize s to 0  
//  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
//    for(int blocknr = 0; blocknr < 4; ++blocknr) {
//      (s[rnd_i][blocknr]).setZero();
//    }
//  }

  // dilution from the left
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
    for(int blocknr = 0; blocknr < 4; ++blocknr) {

      // s holds left diluted basicoperator
      for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
        s[rnd_i][blocknr].row(vec_i % quarks[0].number_of_dilution_E) +=
            std::conj(rnd_vec[rnd_i](blocknr + vec_i * 4 + 
              4 * number_of_eigen_vec * t_sink)) * 
            basicoperator[p][t_sink][displ].row(vec_i);
      }
    }
  }
                
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {

    for(int col = 0; col < 4; ++col) {
      for(int row = 0; row  < 4; ++row){

        // propagator D_u^-1 = perambulator(tsource, tsink) * basicoperator(tsink)
        // calculate columns of D_u^-1. gamma structure can be implented by
        // reordering columns and multiplying them with constants

        contraction[particle_no][p][rnd_i][col].block(row * number_of_eigen_vec, 0,
            number_of_eigen_vec, number_of_eigen_vec) =
          perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + 
            number_of_eigen_vec * row,
            (quarks[0].number_of_dilution_E) * quarks[0].number_of_dilution_D * 
            t_sink_dil + (quarks[0].number_of_dilution_E) * col,
            number_of_eigen_vec,
            (quarks[0].number_of_dilution_E)) *
          s[rnd_i][col];

      }
    }

  }

//  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; rndi++){
//    for(int blocknr = 0; blocknr < 4; blocknr++){
//      delete [] s[rnd_i][blocknr];
//    }
//    delete [] s[rnd_i];
//  }
//  delete [] s;

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  return;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


// initializes contractions_dagger[col] with columns of D_d^-1

void BasicOperator::init_operator_d (const int particle_no, const int t_source, 
    const int t_sink, ReadWrite* rewr, const char dilution, const int p, 
    const int displ){

  clock_t t = clock();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  int t_sink_dil;
  const vec_eigen perambulator = rewr->get_perambulator(); 
  const dim3_eigen_array basicoperator = rewr->get_basicoperator(); 

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

                
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {

    for(int col = 0; col < 4; ++col) {
      for(int row = 0; row  < 4; ++row){

        // propagator D_d^-1 = perambulator(t_source, t_sink)^dagger * 
        // basicoperator(tsource) 
        // = gamma_5 D_u^-1 gamma_5 according to gamma_5 trick
        // only necassary to build this for charged particles.
        // TODO: implement a flag to omit this calculation
        // TODO: implement more versatile momentum structure

        contraction_dagger[particle_no][p][rnd_i][row].block(
            col * quarks[0].number_of_dilution_E, 0,
            quarks[0].number_of_dilution_E, number_of_eigen_vec) = 
          (perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + 
            number_of_eigen_vec * row,
            (quarks[0].number_of_dilution_E) * quarks[0].number_of_dilution_D * 
            t_sink_dil + (quarks[0].number_of_dilution_E) * col,
            number_of_eigen_vec,
            (quarks[0].number_of_dilution_E))).adjoint() *
          basicoperator[p][t_source][displ];
          
        // that's the best criterium I could think up for multiplication with
        // gamma_5 from left and right side. It changes the sign of the two
        // upper right and two lower left blocks in dirac space
        if( ((row + col) == 3) || (abs(row - col) > 1) ){
        contraction_dagger[particle_no][p][rnd_i][row].block(col * 
            quarks[0].number_of_dilution_E, 0,
            quarks[0].number_of_dilution_E, number_of_eigen_vec) *= -1;
        }

      }
    }

  }

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  return;
}

/******************************************************************************/
/******************************************************************************/
/*charged*********************************************************************/
/******************************************************************************/
/******************************************************************************/

// returns D_u^-1 Gamma

void BasicOperator::get_operator_charged (dim2_eigen_array& op_1, 
    const int particle_no, const int t_sink, ReadWrite* rewr, const int dirac, 
    const int p) const {

  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const vec_eigen rnd_vec = rewr->get_random_vector();

  // for charged particles D_u^-1 must be diluted from the right side
  // to match the dilution in the perambulator for D_d^-1

  // s_temp holds diluted contraction matrix
  Eigen::MatrixXcd s_temp = Eigen::MatrixXcd::Zero(number_of_eigen_vec,
      quarks[0].number_of_dilution_E);

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = rnd_i + 1; rnd_j < number_of_rnd_vec; ++rnd_j) {

      // column number of contraction
      for(int col = 0; col < 4; col++){
      
        // row number of contraction corresponds to blocknr in init_operator
        for(int row = 0; row  < 4; row++) {

          // dilutes contraction in blocks of 
          // number_of_eigen_vec * number_of_dilution_E (one dilution step for
          // every dirac and number_of_dilution_D index
          for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
            s_temp.col(vec_i % quarks[0].number_of_dilution_E) +=
              (rnd_vec[rnd_j](col + vec_i * 4 + 
                4 * number_of_eigen_vec * t_sink)) *
              contraction[particle_no][p][rnd_i][col].block(
                row * number_of_eigen_vec, vec_i, number_of_eigen_vec, 1);
          }
  
          // introducing gamma structure via reordering of columns. Blockwise
          // due to different randomvector-entries
          (op_1[rnd_i][rnd_j]).block(row * number_of_eigen_vec, 
              gamma[dirac].row[col] * quarks[0].number_of_dilution_E, 
              number_of_eigen_vec, quarks[0].number_of_dilution_E) = 
            gamma[dirac].value[col] * s_temp;

          // reset s_temp to 0
          (s_temp).setZero();

        } //end for row
  
      } // end for col

    } // end for rnd_j
  } // end for rnd_i

//  delete [] s_temp;
  return;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//returns D_d^-1 Gamma

void BasicOperator::get_operator_g5 (vec_eigen& op_1, 
    const int particle_no, const int dirac, const int p) const{

  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

    for(int i = 0; i < 4; i++) {
      op_1[rnd_i].block(0, gamma[dirac].row[i] * number_of_eigen_vec,
//          4 * number_of_eigen_vec, number_of_eigen_vec) =
          4 * quarks[0].number_of_dilution_E, number_of_eigen_vec) =
      gamma[dirac].value[i] * contraction_dagger[particle_no][p][rnd_i][i];
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

  return;
} 
 
// TODO: think about speedup from extracting factors -1 and +-i in get_operator 
// and get_operator_g5



/******************************************************************************/
/******************************************************************************/
/*uncharged********************************************************************/
/******************************************************************************/
/******************************************************************************/



// returns D_u^-1 Gamma

void BasicOperator::get_operator_uncharged (vec_eigen& op_1, 
    const int particle_no, const int dirac, const int p) const{

  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

    for(int i = 0; i < 4; i++) {
 
      op_1[rnd_i].block(0, gamma[dirac].row[i] * number_of_eigen_vec,
          4 * number_of_eigen_vec, number_of_eigen_vec) =
      gamma[dirac].value[i] * contraction[particle_no][p][rnd_i][i];
  
    }
  }
  return;
}


