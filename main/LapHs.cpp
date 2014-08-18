//============================================================================
// Name        : LapHs.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : stochastic LapH code
//============================================================================

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "GlobalData.h"
#include "BasicOperator.h"
#include "ReadWrite.h"

int main (int ac, char* av[]) {
  // ***************************************************************************
  // ***************************************************************************
  // initialization ************************************************************
  // ***************************************************************************
  // ***************************************************************************

  // initialization of eigen OMP paralization
  Eigen::initParallel();
  // set numer of threads used by eigen
  // first line sets number of threads directly
  // second line lets OMP decide on the number of threads,
  // e. g. via OMP_NUM_THREADS
  //Eigen::setNbThreads(4);
  Eigen::setNbThreads(0);
  //check the number of threads used
  const int nthreads = Eigen::nbThreads();

  std::cout << "contraction code for stochastic dilution" << std::endl;
  std::cout << "using " << nthreads << " threads for eigen\n" << std::endl;

  // reading in global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  // reading in of data
  ReadWrite* rewr = new ReadWrite;

  // everything for operator handling
  BasicOperator* basic = new BasicOperator();

  // global variables from input file needed in main function
  const int Lt = global_data->get_Lt();
  const int end_config = global_data->get_end_config();
  const int delta_config = global_data->get_delta_config();
  const int start_config = global_data->get_start_config();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_squared = number_of_max_mom * number_of_max_mom;
  const int number_of_momenta = global_data->get_number_of_momenta();

  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const std::vector<int> mom_squared = global_data->get_momentum_squared();

  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();

  const int displ_min = global_data->get_displ_min();
  const int displ_max = global_data->get_displ_max();

  const int p_min = 0; //number_of_momenta/2;
  const int p_max = number_of_momenta;

  std::string outpath = global_data->get_output_path() + "/";

  // other variables
  clock_t time;

  const std::complex<double> I(0.0, 1.0);

  char outfile[400];
  FILE *fp = NULL;

  // ***************************************************************************
  // ***************************************************************************
  // memory allocation *********************************************************
  // ***************************************************************************
  // ***************************************************************************

  // memory for operators (old)
//  Eigen::MatrixXcd op_D_tsource = Eigen::MatrixXcd::Zero(
//      quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
//      4 * number_of_eigen_vec);
//  Eigen::MatrixXcd op_D_tsource_1 = Eigen::MatrixXcd::Zero(
//      quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
//      4 * number_of_eigen_vec);
//  Eigen::MatrixXcd op_D_tsink = Eigen::MatrixXcd::Zero(
//      quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
//      4 * number_of_eigen_vec);
//  Eigen::MatrixXcd op_D_tsink_1 = Eigen::MatrixXcd::Zero(
//      quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
//      4 * number_of_eigen_vec);
//  Eigen::MatrixXcd*** X = new Eigen::MatrixXcd**[number_of_rnd_vec];
//  Eigen::MatrixXcd*** Y = new Eigen::MatrixXcd**[number_of_rnd_vec];
//  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
//    X[rnd1] = new Eigen::MatrixXcd*[number_of_rnd_vec];
//    Y[rnd1] = new Eigen::MatrixXcd*[number_of_rnd_vec];
//    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++){
//      X[rnd1][rnd2] = new Eigen::MatrixXcd[number_of_rnd_vec];
//      Y[rnd1][rnd2] = new Eigen::MatrixXcd[number_of_rnd_vec];
//      for(int rnd3 = 0; rnd3 < number_of_rnd_vec; rnd3++){
//        X[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
//            4 * quarks[0].number_of_dilution_E, 4 * quarks[0].number_of_dilution_E);
//        Y[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
//            4 * quarks[0].number_of_dilution_E, 4 * quarks[0].number_of_dilution_E);
//      }
//    }
//  }

  // memory for the correlation function
  std::complex<double>**** C2_mes = 
      new std::complex<double>***[number_of_momenta];
  for(int p_u = 0; p_u < number_of_momenta; ++p_u){
    C2_mes[p_u] = new std::complex<double>**[number_of_momenta];
    for(int p_d = 0; p_d < number_of_momenta; ++p_d){
      C2_mes[p_u][p_d] = new std::complex<double>*[16]; // 16 Dirac matrices
      for(int dirac = 0; dirac < 16; ++dirac){
        C2_mes[p_u][p_d][dirac] = new std::complex<double>[Lt];
      }
    }
  }

//  std::complex<double>**** C2_dis = 
//      new std::complex<double>***[number_of_momenta];
//  for(int p = 0; p < number_of_momenta; ++p){
//    C2_dis[p] = new std::complex<double>**[16]; // 16 Dirac matrices
//    for(int dirac = 0; dirac < 16; ++dirac){
//      C2_dis[p][dirac] = new std::complex<double>*[number_of_rnd_vec];
//      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++) {
//        C2_dis[p][dirac][rnd1] = new std::complex<double>[Lt];
//      }
//    }
//  }
//
//  std::complex<double>***** C4_mes = 
//      new std::complex<double>****[number_of_momenta];
//  for(int p_u = 0; p_u < number_of_momenta; p_u++){
//    C4_mes[p_u] = new std::complex<double>***[number_of_momenta];
//    for(int p_d = 0; p_d < number_of_momenta; p_d++){
//      C4_mes[p_u][p_d] = new std::complex<double>**[16];
//      for(int dirac_u = 0; dirac_u < 16; ++dirac_u){
//        C4_mes[p_u][p_d][dirac_u] = new std::complex<double>*[16];
//        for(int dirac_d = 0; dirac_d < 16; ++dirac_d){
//          C4_mes[p_u][p_d][dirac_u][dirac_d] = new std::complex<double>[Lt];
//        }
//      }
//    }
//  }

  std::complex<double>******** Corr = 
      new std::complex<double>*******[number_of_momenta];
  for(int p1 = 0; p1 < number_of_momenta; ++p1){
    Corr[p1] = new std::complex<double>******[number_of_momenta];
    for(int p2 = 0; p2 < number_of_momenta; ++p2){
      Corr[p1][p2] = new std::complex<double>*****[16];
      for(int dirac1 = 5; dirac1 < 6; ++dirac1){
        Corr[p1][p2][dirac1] = new std::complex<double>****[16];
        for(int dirac2 = 5; dirac2 < 6; ++dirac2){
          Corr[p1][p2][dirac1][dirac2] = new std::complex<double>***[Lt];
          for(int t1 = 0; t1 < Lt; t1++){
            Corr[p1][p2][dirac1][dirac2][t1] = new std::complex<double>**[Lt];
            for(int t2 = 0; t2 < Lt; t2++){
              Corr[p1][p2][dirac1][dirac2][t1][t2] = 
                  new std::complex<double>*[number_of_rnd_vec];
              for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
                Corr[p1][p2][dirac1][dirac2][t1][t2][rnd1] = 
                    new std::complex<double>[number_of_rnd_vec];
              }
            }
          }
        }
      }
    }
  }

  // memory for traces
//  double part1[number_of_rnd_vec][number_of_rnd_vec];
//  double part2[number_of_rnd_vec][number_of_rnd_vec];

//  int norm = 0;
//  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
//    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
//      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
//        if((rnd2 != rnd1) && (rnd2 != rnd3)){
//          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
//            if((rnd4 != rnd1) && (rnd4 != rnd3)){
//              norm++;
//              //std::cout << "\n\nnorm: " << norm << rnd1 << rnd3 << rnd2 << rnd4 << std::endl;
//            }
//          }
//        }
//      }
//    }
//  }

//  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
//  const double norm1 = Lt * norm;

  // memory for operators
  // from t_source to t_sink (op_1) and vice versa (op_2)
  // additional t_source to t_sink (op_3) and t_sink to t_source (op_4) for
  // 4-point functions
  // op_5 and op_6 for neutral particles
  Eigen::MatrixXcd** op_1 = new Eigen::MatrixXcd*[number_of_rnd_vec];
//  Eigen::MatrixXcd** op_3 = new Eigen::MatrixXcd*[number_of_rnd_vec];
  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
    op_1[rnd_i] = new Eigen::MatrixXcd[number_of_rnd_vec];
//    op_3[rnd_i] = new Eigen::MatrixXcd[number_of_rnd_vec];
  }
  Eigen::MatrixXcd* op_2 = new Eigen::MatrixXcd[number_of_rnd_vec];
//  Eigen::MatrixXcd* op_4 = new Eigen::MatrixXcd[number_of_rnd_vec];
//
//  Eigen::MatrixXcd* op_5 = new Eigen::MatrixXcd[number_of_rnd_vec];
//  Eigen::MatrixXcd* op_6 = new Eigen::MatrixXcd[number_of_rnd_vec];

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j){
      op_1[rnd_i][rnd_j] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
          4 * quarks[0].number_of_dilution_E);
//      op_3[rnd_i][rnd_j] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
//          4 * quarks[0].number_of_dilution_E);
    }

    op_2[rnd_i] = Eigen::MatrixXcd(4 * quarks[0].number_of_dilution_E, 
        4 * number_of_eigen_vec);
//    op_4[rnd_i] = Eigen::MatrixXcd(4 * quarks[0].number_of_dilution_E, 
//        4 * number_of_eigen_vec);
//
//    op_5[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
//        4 * number_of_eigen_vec);
//    op_6[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
//        4 * number_of_eigen_vec);

  }

  // ***************************************************************************
  // ***************************************************************************
  // Loop over all configurations **********************************************
  // ***************************************************************************
  // ***************************************************************************

  for(int config_i = start_config; config_i <= end_config; config_i +=
      delta_config){

    std::cout << "\nprocessing configuration: " << config_i << "\n\n";

    rewr->read_perambulators_from_file(config_i);
    rewr->read_rnd_vectors_from_file(config_i);
//    rewr->read_eigenvectors_from_file(config_i);
//    rewr->read_lime_gauge_field_doubleprec_timeslices(config_i);
    rewr->build_source_matrix(config_i);


    // *************************************************************************
    // TWO PT CONTRACTION 1 ****************************************************
    // *************************************************************************

    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of pi_+/-:\n";
    time = clock();

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac = 0; dirac < 16; ++dirac)
          for(int t1 = 0; t1 < Lt; ++t1)
            C2_mes[p_u][p_d][dirac][t1] = std::complex<double>(0.0, 0.0);

    for(int p1 = 0; p1 < number_of_momenta; ++p1)
      for(int p2 = 0; p2 < number_of_momenta; ++p2)
        for(int dirac1 = 5; dirac1 < 6; ++dirac1)
          for(int dirac2 = 5; dirac2 < 6; ++dirac2)
            for(int t1 = 0; t1 < Lt; ++t1)
              for(int t2 = 0; t2 < Lt; ++t2)
                for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++)
                  for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++)
                    Corr[p1][p2][dirac1][dirac2][t1][t2][rnd1][rnd2] = 
                        std::complex<double>(0.0, 0.0);


//    for(int p = 0; p < number_of_momenta; ++p)
//      for(int dirac = 0; dirac < 16; ++dirac)
//        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++)
//          for(int t1 = 0; t1 < Lt; ++t1)
//            C2_dis[p][dirac][rnd1][t1] = std::complex<double>(0.0, 0.0);

#if 0
    // disconnected part of pi^0
    for(int t_source = 0; t_source < Lt; ++t_source){

      for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
        for(int p = p_min; p < p_max; ++p) {

          basic->init_operator_uncharged(t_source, t_source, rewr, 'b', p, 0);
          basic->get_operator_uncharged(op_5, dirac);
    
          for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
             C2_dis[p][dirac][rnd1][t_source] += 
                 (op_5[rnd1]).trace();
          }

        }
      }

    }

#endif
#if 1
    // pi^+-
    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        for(int p = p_min; p < p_max; ++p) {
          // initialize contraction[rnd_i] = perambulator * basicoperator
          // = D_u^-1
          // choose 'i' for interlace or 'b' for block time dilution scheme
          basic->init_operator_u(0, t_source, t_sink, rewr, 'b', p, 0);
          basic->init_operator_d(0, t_source, t_sink, rewr, 'b', p, 0);
        }

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
          for(int p_u = p_min; p_u < p_max; ++p_u) {

            // "multiply contraction[rnd_i] with gamma structure"
            // contraction[rnd_i] are the columns of D_u^-1 which get
            // reordered by gamma multiplication. No actual multiplication
            // is carried out
            basic->get_operator_charged(op_1, 0, t_sink, rewr, dirac_u, p_u);

            for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
              for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){
  
                  // same as get_operator but with gamma_5 trick. D_u^-1 is
                  // daggered and multipied with gamma_5 from left and right
                  // the momentum is changed to reflect the switched sign in
                  // the momentum exponential for pi_+-
                  basic->get_operator_g5(op_2, 0, dirac_d, p_d);
     
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
                      // building Correlation function get quite intuitive
                      // C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                      // TODO: find signflip of imaginary part
                      // TODO: is C2_mes[dirac][p] better?
                      Corr[p_u][p_d][dirac_u][dirac_d][t_source][t_sink][rnd1][rnd2] = 
                        (op_2[rnd2] * op_1[rnd1][rnd2]).trace();
    
    
                    } // end for rnd2
                  } // end for rnd1 
  
                } // end if momentum squared
              } // end for momentum d-quark
            } // end for dirac d-quark

          } // end for momentum u-quark
        } // end for dirac u-quark

      } // end for t_sink
    } // end for t_source

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
          for(int p_u = p_min; p_u < p_max; ++p_u) {
            for(int p_d = p_min; p_d < p_max; ++p_d) {
    
              for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
      
                  C2_mes[p_u][p_d][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
                      Corr[p_u][number_of_momenta - p_d - 1][dirac][dirac]
                      [t_source][t_sink][rnd1][rnd2];
                }
              }
      
            }
          }
        }

      }
    }
#endif
#if 0
// code for pi0 

            basic->init_operator_uncharged(t_source, t_sink, rewr, 'b', p, 0);
            basic->get_operator_uncharged(op_5, dirac);
    
            basic->init_operator_uncharged(t_sink, t_source, rewr, 'b', 
                number_of_momenta - p - 1, 0);
            basic->get_operator_uncharged(op_6, dirac);            

            for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
              for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
                  // building Correlation function get quite intuitive
                  // C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                  // TODO: find signflip of imaginary part
                  // TODO: is C2_mes[dirac][p] better?
                  C2_mes[p][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
                      (op_6[rnd2] * op_5[rnd1]).trace();

              }
            }   

          }
        }

      }
    }
#endif
    // normalization of correlation function
    double norm3 = Lt * number_of_rnd_vec * (number_of_rnd_vec - 1) * 0.5;
    for(int t = 0; t < Lt; ++t){
      for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
        for(int p_u = p_min; p_u < p_max; ++p_u) {
          for(int p_d = p_min; p_d < p_max; ++p_d) {
            C2_mes[p_u][p_d][dirac][t] /= norm3;
          }
        }
      }
    }

    // output to binary file
//    sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_pi+-_conf%04d.dat", 
//        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, 
//        displ_min, displ_max, config_i);
    sprintf(outfile, "%s/C2_pi+-_conf%04d.dat", outpath.c_str(), config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    else {
      for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
        for(int p = 0; p <= max_mom_squared; p++){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){
              fwrite((double*) C2_mes[p_u][p_u][dirac], sizeof(double), 2 * Lt, fp);
            }
          }
        }
#if 0
        for(int p = 1; p <= max_mom_squared; p++){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){
              fwrite((double*) C2_mes[number_of_momenta / 2][p_u][dirac], sizeof(double), 2 * Lt, fp);
            }
          }
        }
#endif
      }
      fclose(fp);
    }

#if 0
    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_pi+-_conf%04d.dat", 
        outpath.c_str(), dirac_min, dirac_max, 0, 
        displ_min, displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
      fwrite((double*) C2_mes[number_of_momenta/2]
          [number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
    fclose(fp);

    for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
      sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_dis_u_rnd%02d_conf%04d.dat", 
          outpath.c_str(), dirac_min, dirac_max, number_of_max_mom, displ_min, 
          displ_max, rnd_i, config_i);
      if((fp = fopen(outfile, "wb")) == NULL)
        std::cout << "fail to open outputfile" << std::endl;
      for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
        for(int p = 0; p < number_of_momenta; ++p)
          fwrite((double*) C2_dis[p][dirac][rnd_i], sizeof(double), 2 * Lt, fp);
      fclose(fp);
    }
#endif

#if 0
    // output to terminal
    printf("\n");
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
      printf("\tdirac    = %02d\n", dirac);
      for(int p = 0; p <= max_mom_squared; p++){
        printf("\tmomentum_u = %02d\n", p);
        printf("\tmomentum_d = %02d\n", p);
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if((mom_squared[p_u] == p)){
            //printf(
            //    "\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
//                  for(int t1 = 0; t1 < Lt; ++t1){
//                    printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
//                        imag(C2_mes[p_u][p_u][dirac][t1]));
//                  }
            printf("\n");
            printf("p_u = %02d\n", p_u);
          }
        }
      }

      for(int p = 1; p <= max_mom_squared; p++){
        printf("\tmomentum_u = %02d\n", 0);
        printf("\tmomentum_d = %02d\n", p);
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if((mom_squared[p_u] == p)){
            //printf(
            //    "\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
//                  for(int t1 = 0; t1 < Lt; ++t1){
//                    printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
//                        imag(C2_mes[p_u][p_u][dirac][t1]));
//                  }
            printf("\n");
            printf("p_u = %02d\n", p_u);
          }
        }
      }

    }
#endif
    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);


    // *************************************************************************
    // FOUR PT CONTRACTION 1 ***************************************************
    // *************************************************************************
#if 0
    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
    time = clock();

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < 16; ++dirac_u)
          for(int dirac_d = 0; dirac_d < 16; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
          for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
             for(int p_u = p_min; p_u < p_max; ++p_u) {
               for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){

  
      //            basic->init_operator_charged(t_source_1, t_sink_1, rewr, 'b', p, 0);
      //    
      //            basic->get_operator_charged(op_1, rewr, dirac, t_sink_1);
      //            basic->get_operator_g5(op_2, dirac);
      //    
      //            basic->init_operator_charged(t_source, t_sink, rewr, 'b', 
      //                number_of_momenta - p - 1, 0);
      //    
      //            basic->get_operator_charged(op_3, rewr, dirac, t_sink);
      //            basic->get_operator_g5(op_4, dirac);
      //    
      //            // first trace
      //            for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
      //              // first u quark: t_source_1 -> t_sink_1
      //              for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
      //                // first d quark: t_sink_1 -> t_source_1
      //                part1[rnd1][rnd3] = std::real((op_2[rnd3] * op_1[rnd1][rnd3]).trace());
      //              }
      //            }
      //    
      //            // second trace
      //            for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
      //              // second u quark: t_source -> t_sink
      //              for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
      //                // second d quark: t_sink -> t_source
      //                part2[rnd2][rnd4] = std::real((op_4[rnd4] * op_3[rnd2][rnd4]).trace());
      //              }
      //            }
          
                  // complete diagramm
                  // every quark line must have its own random vec
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
                      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
                        if((rnd2 != rnd1) && (rnd2 != rnd3)){
                          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                            if((rnd4 != rnd1) && (rnd4 != rnd3)){
                              C4_mes[p_u][p_d][dirac_u][dirac_d]
                                  [abs((t_sink - t_source - Lt) % Lt)] +=
                                (Corr[p_u]
                                  [number_of_momenta - p_d - 1]
                                  [dirac_u][dirac_d]
                                  [t_source_1][t_sink_1][rnd1][rnd3]) *
                                (Corr[number_of_momenta - p_u - 1]
                                  [p_d][dirac_u][dirac_d]
                                  [t_source][t_sink][rnd2][rnd4]);
      //                        C2_mes[p][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
      //                            part1[rnd1][rnd3] * part2[rnd2][rnd4];
                            }
                          }
                        }
                      }
                    }
                  }

                }
              }
            }    
          }
        }

      }
    }

    for(int t = 0; t < Lt; ++t){
      for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
        for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
           for(int p_u = p_min; p_u < p_max; ++p_u) {
             for(int p_d = p_min; p_d < p_max; ++p_d) {
              if(mom_squared[p_u] <= mom_squared[p_d]){
                C4_mes[p_u][p_d][dirac_u][dirac_d][t] /= norm1;
              }
            }
          }
        }
      }
    }


    // output to binary file
    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_1_conf%04d_diag.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int p = 0; p <= max_mom_squared; p++){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){
              fwrite((double*) C4_mes[p_u][p_u][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
            }
          }
        }
      }
    }
    fclose(fp);

    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_1_conf%04d_all.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int offset = 0; offset <= max_mom_squared; offset++){
          for(int p = 0; p <= max_mom_squared; p++){
            if((p + offset) <= max_mom_squared){
              for(int p_u = p_min; p_u < p_max; ++p_u){
                if(mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(mom_squared[p_d] == (p + offset)){
                      fwrite((double*) C4_mes[p_u][p_d][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    fclose(fp);

//    sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_1_conf%04d.dat", 
//        outpath.c_str(), dirac_min, dirac_max, 0, displ_min, 
//        displ_max, config_i);
//    if((fp = fopen(outfile, "wb")) == NULL)
//      std::cout << "fail to open outputfile" << std::endl;
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
//      fwrite((double*) C2_mes[number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
//    fclose(fp);

    // output to terminal
    printf("\n");
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
      printf("\tdirac    = %02d\n", dirac);
      for(int offset = 0; offset <= max_mom_squared; offset++){
        for(int p = 0; p <= max_mom_squared; p++){
          if((p + offset) <= max_mom_squared){
            printf("\tmomentum_u = %02d\n", p);
            printf("\tmomentum_d = %02d\n", p + offset);
            for(int p_u = p_min; p_u < p_max; ++p_u){
              if((mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(mom_squared[p_d] == (p + offset)){
                    //printf(
                    //    "\t t\tRe(C4_1_con)\tIm(C4_1_con)\n\t----------------------------------\n");
//                    for(int t1 = 0; t1 < Lt; ++t1){
//                      printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p_u][p_d][dirac][dirac][t1]),
//                          imag(C4_mes[p_u][p_d][dirac][dirac][t1]));
//                    }
                  printf("\n");
                  printf("p_u = %02d\tp_d = %02d\n", p_u, p_d);
                  }
                }
              }
            }
          }
        printf("\n");
        }
      }
    printf("\n");
    }

    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
#endif
    // *************************************************************************
    // FOUR PT CONTRACTION 2 ***************************************************
    // *************************************************************************
#if 0
    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
    time = clock();

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < 16; ++dirac_u)
          for(int dirac_d = 0; dirac_d < 16; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt - 1; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
          for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
             for(int p_u = p_min; p_u < p_max; ++p_u) {
               for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){


      //            basic->init_operator_charged(t_source_1, t_sink, rewr, 'b', p, 0);
      //    
      //            basic->get_operator_charged(op_1, rewr, dirac, t_sink);
      //            basic->get_operator_g5(op_2, dirac);
      //    
      //            basic->init_operator_charged(t_source, t_sink_1, rewr, 'b', 
      //                number_of_momenta - p - 1, 0);
      //    
      //            basic->get_operator_charged(op_3, rewr, dirac, t_sink_1);
      //            basic->get_operator_g5(op_4, dirac);
      //    
      //            // first trace
      //            for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
      //              // first u quark: t_source_1 -> t_sink
      //              for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
      //                // first d quark: t_sink -> t_source_1
      //                part1[rnd1][rnd3] = std::real((op_2[rnd3] * op_1[rnd1][rnd3]).trace());
      //              }
      //            }
      //    
      //            // second trace
      //            for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
      //              // second u quark: t_source -> t_sink_1
      //              for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
      //                // second d quark: t_sink_1 -> t_source
      //                part2[rnd2][rnd4] = std::real((op_4[rnd4] * op_3[rnd2][rnd4]).trace());
      //              }
      //            }
          
                  // complete diagramm
                  // every quark line must have its own random vec
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
                      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
                        if((rnd2 != rnd1) && (rnd2 != rnd3)){
                          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                            if((rnd4 != rnd1) && (rnd4 != rnd3)){
                              C4_mes[p_u][p_d][dirac_u][dirac_d]
                                  [abs((t_sink - t_source - Lt) % Lt)] +=
                                (Corr[p_u][number_of_momenta - p_d - 1]
                                  [dirac_u][dirac_d][t_source_1][t_sink][rnd1][rnd3]) *
                                (Corr[number_of_momenta - p_u - 1][p_d]
                                  [dirac_u][dirac_d][t_source][t_sink_1][rnd2][rnd4]);
                            }
                          }
                        }
                      }
                    }
                  }

                }
              }
            }
          }
        }

      }
    }

    for(int t = 0; t < Lt; ++t){
      for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
        for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
           for(int p_u = p_min; p_u < p_max; ++p_u) {
             for(int p_d = p_min; p_d < p_max; ++p_d) {
              if(mom_squared[p_u] <= mom_squared[p_d]){
                C4_mes[p_u][p_d][dirac_u][dirac_d][t] /= norm1;
              }
            }
          }
        }
      }
    }

    // output to binary file
    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_2_conf%04d_diag.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
          std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int p = 0; p <= max_mom_squared; p++){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){
              fwrite((double*) C4_mes[p_u][p_u][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
            }
          }
        }
      }
    }
    fclose(fp);

    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_2_conf%04d_all.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
          std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int offset = 0; offset <= max_mom_squared; offset++){
          for(int p = 0; p <= max_mom_squared; p++){
            if((p + offset) <= max_mom_squared){
              for(int p_u = p_min; p_u < p_max; ++p_u){
                if(mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(mom_squared[p_d] == (p + offset)){
                      fwrite((double*) C4_mes[p_u][p_d][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    fclose(fp);

//    sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_2_conf%04d.dat", 
//        outpath.c_str(), dirac_min, dirac_max, 0, displ_min, 
//        displ_max, config_i);
//    if((fp = fopen(outfile, "wb")) == NULL)
//      std::cout << "fail to open outputfile" << std::endl;
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
//      fwrite((double*) C2_mes[number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
//    fclose(fp);

    // output to terminal
    printf("\n");
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
      printf("\tdirac    = %02d\n", dirac);
      for(int offset = 0; offset <= max_mom_squared; offset++){
        for(int p = 0; p <= max_mom_squared; p++){
          if((p + offset) <= max_mom_squared){
            printf("\tmomentum_u = %02d\n", p);
            printf("\tmomentum_d = %02d\n", p + offset);
            for(int p_u = p_min; p_u < p_max; ++p_u){
              if((mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(mom_squared[p_d] == (p + offset)){
                    //printf(
                    //    "\t t\tRe(C4_2_con)\tIm(C4_2_con)\n\t----------------------------------\n");
//                    for(int t1 = 0; t1 < Lt; ++t1){
//                      printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
//                          imag(C4_mes[p][p][dirac][dirac][t1]));
//                    }
                    printf("\n");
                    printf("p_u = %02d\tp_d = %02d\n", p_u, p_d);
                  }
                }
              }
            }
          }
          printf("\n");
        }
      }
      printf("\n");
    }

    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
#endif

    // *************************************************************************
    // FOUR PT CONTRACTION 3 ***************************************************
    // *************************************************************************
#if 0
    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of C4_3:\n";
    time = clock();
    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < 16; ++dirac_u)
          for(int dirac_d = 0; dirac_d < 16; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        std::cout << "\tt_source = " << t_source << "\tt_sink = " << t_sink << std::endl;


        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            basic->init_operator_u(0, t_source, t_sink, rewr, 'b', p_u, 0);
            basic->init_operator_u(1, t_source_1, t_sink_1, rewr, 'b', 
                number_of_momenta - p_u - 1, 0);
            }
        }


        for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
          for(int p_d = p_min; p_d < p_max; ++p_d){
            basic->init_operator_d(0, t_source_1, t_sink, rewr, 'b', p_d, 0);
            basic->init_operator_d(1, t_source, t_sink_1, rewr, 'b', 
                number_of_momenta - p_d - 1, 0);
          }
        }

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
           for(int p_u = p_min; p_u < p_max; ++p_u) {
            for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
               for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){

                  basic->get_operator_g5(op_2, 0, dirac_d, 
                      number_of_momenta - p_d - 1);
                  basic->get_operator_charged(op_3, 1, t_sink_1, rewr, dirac_u, 
                      number_of_momenta - p_u - 1);
      
                  // first part
                  for(int rnd3 = 1; rnd3 < number_of_rnd_vec; ++rnd3){
                    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
                      if(rnd2 != rnd3){
                        // second u quark: t_source_1 -> t_sink_1
                        for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                          if(rnd4 != rnd3){
                            // first d quark: t_sink_1 -> t_source
                            X[rnd3][rnd2][rnd4] = op_2[rnd3] * op_3[rnd2][rnd4] ;
                          }
                        }
                      }
                    }
                  }

                  basic->get_operator_g5(op_4, 1, dirac_d, p_d);
                  basic->get_operator_charged(op_1, 0, t_sink, rewr, dirac_u, p_u);

                  // second part
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){      
                      // first u quark: t_source -> t_sink
                      for(int rnd4 = 1; rnd4 < number_of_rnd_vec; ++rnd4){
                        if((rnd4 != rnd1) && (rnd4 != rnd3)){
                          // second d quark: t_sink -> t_source_1
                          Y[rnd4][rnd1][rnd3] = op_4[rnd4] * op_1[rnd1][rnd3] ;
                        }
                      }
                    }
                  }
          
                  // complete diagramm
                  // every quark line must have its own random vec
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
                      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
                        if((rnd2 != rnd1) && (rnd2 != rnd3)){
                          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                            if((rnd4 != rnd1) && (rnd4 != rnd3)){
//                              std::cout << "rnd 4 1 3\t" << rnd4 << rnd1 << rnd3 << std::endl;
//                              std::cout << "\n" << Y[rnd4][rnd1][rnd3] << "\n" << std::endl;
//                              std::cout << "rnd 3 2 4\t" << rnd3 << rnd2 << rnd4 << std::endl;
//                              std::cout << "\n" << X[rnd3][rnd2][rnd4] << "\n" << std::endl;
                              C4_mes[p_u][p_d][dirac_u][dirac_d][abs((t_sink - t_source - Lt) % Lt)] += 
                                  std::real((X[rnd3][rnd2][rnd4] * Y[rnd4][rnd1][rnd3]).trace());
                            }
                          }
                        }
                      }
                    }
                  }

                }
              }
            }
          }
        }

      }
    }

    for(int t = 0; t < Lt; ++t)
      for(int p = 0; p < number_of_momenta; ++p)
        for(int p1 = 0; p1 < number_of_momenta; ++p1)
        for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
          C4_mes[p][p1][dirac][dirac][t] /= norm1;

    // output to binary file
    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_3_conf%04d_diag.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int p = 0; p <= max_mom_squared; p++){
          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){
              fwrite((double*) C4_mes[p_u][p_u][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
            }
          }
        }
      }
    }
    fclose(fp);

    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_3_conf%04d_all.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
        for(int offset = 0; offset <= max_mom_squared; offset++){
          for(int p = 0; p <= max_mom_squared; p++){
            if((p + offset) <= max_mom_squared){
              for(int p_u = p_min; p_u < p_max; ++p_u){
                if(mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(mom_squared[p_d] == (p + offset)){
                      fwrite((double*) C4_mes[p_u][p_d][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    fclose(fp);

#if 0
    sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_3_conf%04d.dat", 
        outpath.c_str(), dirac_min, dirac_max, 0, displ_min, 
        displ_max, config_i);
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile" << std::endl;
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
      fwrite((double*) C4_mes[number_of_momenta/2]
          [number_of_momenta/2][dirac][dirac], sizeof(double), 2 * Lt, fp);
    fclose(fp);
#endif

    // output to terminal
    printf("\n");
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
      printf("\tdirac    = %02d\n", dirac);
      for(int p = p_min; p < p_max; ++p) {
        printf("\tmomentum = %02d\n", p);
        //printf(
        //    "\t t\tRe(C4_3_con)\tIm(C4_3_con)\n\t----------------------------------\n");
        for(int t1 = 0; t1 < Lt; ++t1){
          printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
              imag(C4_mes[p][p][dirac][dirac][t1]));
        }
        printf("\n");
      }
      printf("\n");
    }


    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
#endif
    // *************************************************************************
    // FOUR PT CONTRACTION 4 ***************************************************
    // *************************************************************************

    // identical to FOUR PT CONTRACTION 3

  } // loop over configs ends here

  // TODO: freeing all memory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for(int p_u = 0; p_u < number_of_momenta; ++p_u){
    for(int p_d = 0; p_d < number_of_momenta; ++p_d){
      for(int dirac = 0; dirac < 16; ++dirac){
        delete [] C2_mes[p_u][p_d][dirac];
      }
      delete [] C2_mes[p_u][p_d];
    }
    delete [] C2_mes[p_u];
  }
  delete [] C2_mes;

//  for(int p = 0; p < number_of_momenta; ++p){
//    for(int dirac = 0; dirac < 16; ++dirac){
//      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++) {
//        delete [] C2_dis[p][dirac][rnd1];
//      }
//      delete [] C2_dis[p][dirac];
//    }
//    delete [] C2_dis[p];
//  }
//  delete [] C2_dis;

//  for(int p_u = 0; p_u < number_of_momenta; p_u++){
//    for(int p_d = 0; p_d < number_of_momenta; p_d++){
//      for(int dirac_u = 0; dirac_u < 16; ++dirac_u){
//        for(int dirac_d = 0; dirac_d < 16; ++dirac_d){
//          delete [] C4_mes[p_u][p_d][dirac_u][dirac_d];
//        }
//        delete [] C4_mes[p_u][p_d][dirac_u];
//      }
//      delete [] C4_mes[p_u][p_d];
//    }
//    delete [] C4_mes[p_u];
//  }
//  delete [] C4_mes;

  for(int p1 = 0; p1 < number_of_momenta; ++p1){
    for(int p2 = 0; p2 < number_of_momenta; ++p2){
      for(int dirac1 = 5; dirac1 < 6; ++dirac1){
        for(int dirac2 = 5; dirac2 < 6; ++dirac2){
          for(int t1 = 0; t1 < Lt; t1++){
            for(int t2 = 0; t2 < Lt; t2++){
              for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
                delete [] Corr[p1][p2][dirac1][dirac2][t1][t2][rnd1];
              }
              delete [] Corr[p1][p2][dirac1][dirac2][t1][t2];
            }
            delete [] Corr[p1][p2][dirac1][dirac2][t1];
          }
          delete [] Corr[p1][p2][dirac1][dirac2];
        }
        delete [] Corr[p1][p2][dirac1];
      }
      delete [] Corr[p1][p2];
    }
    delete [] Corr[p1];
  }
  delete [] Corr;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
    delete [] op_1[rnd_i];
//    delete [] op_3[rnd_i];
  }
  delete [] op_1;
  delete [] op_2;
//  delete [] op_3;
//  delete [] op_4;
//  delete [] op_5;
//  delete [] op_6;

  delete rewr;
  delete basic;
}

