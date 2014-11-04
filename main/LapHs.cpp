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

#include "GlobalData.h"
#include "BasicOperator.h"
#include "ReadWrite.h"
#include "typedefs.h"

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
  ReadWrite rewr;
  LapH::VdaggerV vdaggerv;

  // everything for operator handling
  BasicOperator* basic = new BasicOperator();

  // global variables from input file needed in main function
  const int Lt = global_data->get_Lt();
  const int end_config = global_data->get_end_config();
  const int delta_config = global_data->get_delta_config();
  const int start_config = global_data->get_start_config();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();

  //const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_squared = global_data->get_number_of_max_mom();
    //number_of_max_mom * number_of_max_mom;
  const int number_of_momenta = global_data->get_number_of_momenta();
  const std::vector<int> mom_squared = global_data->get_momentum_squared();

  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();
  std::vector<int> dirac_ind {5};
  const int number_of_dirac = dirac_ind.size();

  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const int number_of_displ = displ_max - displ_min + 1;

  const int p_min = 0; //number_of_momenta/2;
  const int p_max = number_of_momenta;

  // TODO decide on path
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

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

  // abbreviations for clearer memory allocation. Wont be used in loops and 
  // when building the contractions

  const size_t nmom = number_of_momenta;
  const size_t nrnd = number_of_rnd_vec;
  const size_t ndir = number_of_dirac;
  const size_t ndis = number_of_displ;

  // memory for intermediate matrices when building C4_3 (save multiplications)

  array_Xcd_d3_eigen X(boost::extents[nrnd][nrnd][nrnd]);
  array_Xcd_d3_eigen Y(boost::extents[nrnd][nrnd][nrnd]);
	for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
		for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++){
      for(int rnd3 = 0; rnd3 < number_of_rnd_vec; rnd3++){
  			X[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
  					4 * quarks[0].number_of_dilution_E, 
            4 * quarks[0].number_of_dilution_E);
  			Y[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
  					4 * quarks[0].number_of_dilution_E, 
            4 * quarks[0].number_of_dilution_E);
      }
		}
	}

  // memory for the correlation function
  array_cd_d7 C2_mes(boost::extents[nmom][nmom][ndir][ndir][ndis][ndis][Lt]);
  //TODO: dont need the memory for p_u^2 > p_d^2
  array_cd_d10 Corr(boost::extents[nmom][nmom][ndir][ndir][ndis][ndis][Lt][Lt][nrnd][nrnd]);

  array_cd_d5 C4_mes(boost::extents[nmom][nmom][ndir][ndir][Lt]);

  int norm = 0;
  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
    for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
        if((rnd2 != rnd1) && (rnd2 != rnd3)){
          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
            if((rnd4 != rnd1) && (rnd4 != rnd3)){
              norm++;
              //std::cout << "\n\nnorm: " << norm << rnd1 << rnd3 << rnd2 << rnd4 << std::endl;
            }
          }
        }
      }
    }
  }

  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  // Memory for propagation matrices (is that a word?) from t_source to t_sink
  // (op_1) and vice versa (op_2)
  // additional t_source to t_sink (op_3) and t_sink to t_source (op_4) for
  // 4-point functions
  // 1, 3 -> u-quarks; 2, 4 -> d-quarks; 5, 6 -> u quarks for neutral particle

  array_Xcd_d2_eigen op_1(boost::extents[nrnd][nrnd]);
  array_Xcd_d2_eigen op_3(boost::extents[nrnd][nrnd]);
  vec_Xcd_eigen op_2(number_of_rnd_vec);
  vec_Xcd_eigen op_4(number_of_rnd_vec);
  vec_Xcd_eigen op_5(number_of_rnd_vec);
  vec_Xcd_eigen op_6(number_of_rnd_vec);

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j){
      op_1[rnd_i][rnd_j] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
          4 * quarks[0].number_of_dilution_E);
      op_3[rnd_i][rnd_j] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
          4 * quarks[0].number_of_dilution_E);
    }

    op_2[rnd_i] = Eigen::MatrixXcd(4 * quarks[0].number_of_dilution_E, 
        4 * number_of_eigen_vec);
    op_4[rnd_i] = Eigen::MatrixXcd(4 * quarks[0].number_of_dilution_E, 
        4 * number_of_eigen_vec);

//    op_5[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
//        4 * number_of_eigen_vec);
//    op_6[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
//        4 * number_of_eigen_vec);

  }

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    for(int rnd_j = 0; rnd_j < number_of_rnd_vec; ++rnd_j){

      op_1[rnd_i][rnd_j] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
          4 * quarks[0].number_of_dilution_E);
       }

    op_2[rnd_i] = Eigen::MatrixXcd(4 * quarks[0].number_of_dilution_E, 
        4 * number_of_eigen_vec);
    
  }

  // ***************************************************************************
  // ***************************************************************************
  // Loop over all configurations **********************************************
  // ***************************************************************************
  // ***************************************************************************

  for(int config_i = start_config; config_i <= end_config; config_i +=
      delta_config){

    std::cout << "\nprocessing configuration: " << config_i << "\n\n";

    rewr.read_perambulators_from_file(config_i);
    rewr.read_rnd_vectors_from_file(config_i);
    vdaggerv.build_source_matrix(config_i);


    // *************************************************************************
    // TWO PT CONTRACTION 1 ****************************************************
    // *************************************************************************

    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of pi_+/-:\r";
    time = clock();

    // setting the correlation function to zero
    for(int p1 = 0; p1 < number_of_momenta; ++p1)
      for(int p2 = 0; p2 < number_of_momenta; ++p2)
        for(int dirac1 = 0; dirac1 < number_of_dirac; ++dirac1)
          for(int dirac2 = 0; dirac2 < number_of_dirac; ++dirac2)
            for(int displ1 = 0; displ1 < number_of_displ; ++displ1)
              for(int displ2 = 0; displ2 < number_of_displ; ++displ2)
                for(int t1 = 0; t1 < Lt; ++t1)
                   for(int t1 = 0; t1 < Lt; ++t1)
                     C2_mes[p1][p2][dirac1][dirac2][displ1][displ2][t1] = std::complex<double>(0.0, 0.0);

    for(int p1 = 0; p1 < number_of_momenta; ++p1)
      for(int p2 = 0; p2 < number_of_momenta; ++p2)
        for(int dirac1 = 0; dirac1 < number_of_dirac; ++dirac1)
          for(int dirac2 = 0; dirac2 < number_of_dirac; ++dirac2)
            for(int displ1 = 0; displ1 < number_of_displ; ++displ1)
              for(int displ2 = 0; displ2 < number_of_displ; ++displ2)
                for(int t1 = 0; t1 < Lt; ++t1)
                  for(int t2 = 0; t2 < Lt; ++t2)
                    for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++)
                      for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++)
                        Corr[p1][p2][dirac1][dirac2][displ1][displ2][t1][t2][rnd1][rnd2] = 
                            std::complex<double>(0.0, 0.0);

    // pi^+/-
    // initializing of Corr: calculate all two-operator traces of the form tr(u \Gamma \bar{d})
    // build all combinations of momenta, dirac_structures and displacements as specified in
    // infile

    for(int displ_u = 0; displ_u < number_of_displ; displ_u++){
      for(int displ_d = 0; displ_d < number_of_displ; displ_d++){

          for(int t_source = 0; t_source < Lt; ++t_source){

            std::cout << "\tcomputing the connected contribution of pi_+/-: " 
                << std::setprecision(2) << (float) t_source/Lt*100 <<"%\r" << std::flush;

            for(int t_sink = 0; t_sink < Lt; ++t_sink){
  
              for(int p = p_min; p < p_max; ++p) {
                // initialize contraction[rnd_i] = perambulator * basicoperator
                // = D_u^-1
                // choose 'i' for interlace or 'b' for block time dilution scheme
                // TODO: get that from input file
                // choose 'c' for charged or 'u' for uncharged particles
                basic->init_operator_u(0, t_source, t_sink, rewr, vdaggerv, 'b', p, displ_min + displ_u);
                basic->init_operator_d(0, t_source, t_sink, rewr, vdaggerv, 'b', p, displ_min + displ_d);
              }
    
              for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
                for(int p_u = p_min; p_u < p_max; ++p_u) {
                  // code for pi+-
    
                  // "multiply contraction[rnd_i] with gamma structure"
                  // contraction[rnd_i] are the columns of D_u^-1 which get
                  // reordered by gamma multiplication. No actual multiplication
                  // is carried out
                  basic->get_operator_charged(op_1, 0, t_sink, rewr, dirac_ind.at(dirac_u), p_u);
    
                  for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
                    for(int p_d = p_min; p_d < p_max; ++p_d) {
                      if(mom_squared[p_u] <= mom_squared[p_d]){
      
                        // same as get_operator but with gamma_5 trick. D_u^-1 is
                        // daggered and multipied with gamma_5 from left and right
                        // the momentum is changed to reflect the switched sign in
                        // the momentum exponential for pi_+-
                        basic->get_operator_g5(op_2, 0, dirac_ind.at(dirac_d), p_d);
           
                        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                          for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){

                            // build all 2pt traces leading to C2_mes
                            // Corr = tr(D_d^-1(t_sink) Gamma 
                            //     D_u^-1(t_source) Gamma)

                            Corr[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d]
                                [t_source][t_sink][rnd1][rnd2] = 
                              (op_2[rnd2] * op_1[rnd1][rnd2]).trace();

//                            std::cout << "p" << p_u << p_d << "dirac" << dirac_u << dirac_d << "\nCorr "
//                                << Corr[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d]
//                                [t_source][t_sink][rnd1][rnd2] << std::endl;
          
          
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

    // build 2pt-function C2_mes for pi^+ from Corr. Equivalent two just summing
    // up traces with same time difference between source and sink (all to all)
    // for every dirac structure, momentum, displacement

    // build 2pt-function C2_mes for pi^+ from Corr. Equivalent two just summing
    // up traces with same time difference between source and sink (all to all)
    // for every dirac structure, momentum, displacement

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

         for(int p_u = p_min; p_u < p_max; ++p_u) {
           for(int p_d = p_min; p_d < p_max; ++p_d) {
            if(mom_squared[p_u] <= mom_squared[p_d]){
              for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
//                int dirac_d = dirac_u;
                for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
                  for(int displ_u = 0; displ_u < number_of_displ; displ_u++){
                    for(int displ_d = 0; displ_d < number_of_displ; displ_d++){
            
                      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                        for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
              
                              // building Correlation function 
                              // C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                              // TODO: find signflip of imaginary part
                              // TODO: is C2_mes[dirac][p] better?

                          C2_mes[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d]
                              [abs((t_sink - t_source - Lt) % Lt)] += 
                            Corr[p_u][number_of_momenta - p_d - 1]
                              [dirac_u][dirac_d][displ_u][displ_d]
                              [t_source][t_sink][rnd1][rnd2];
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

    // normalization of correlation function
    double norm3 = Lt * number_of_rnd_vec * (number_of_rnd_vec - 1) * 0.5;
    for(int p_u = p_min; p_u < p_max; ++p_u) {
      for(int p_d = p_min; p_d < p_max; ++p_d) {
        if(mom_squared[p_u] <= mom_squared[p_d]){
          for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
            for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
              for(int displ_u = 0; displ_u < number_of_displ; ++displ_u){
                for(int displ_d = 0; displ_d < number_of_displ; ++displ_d){
                  for(int t = 0; t < Lt; ++t){
                    C2_mes[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d][t] /= norm3;
                  }
                }
              }
            }
          }
        }
      }
    }


    std::cout << "\tcomputing the connected contribution of pi_+/-: 100.00%\n" << std::endl;

    // output to binary file
    // to build a GEVP, the correlators are written into a seperate folder
    // for every dirac structure, momentum, displacement (entry of the GEVP
    // matrix). In the folders a file is created for every configuration which
    // contains all momentum combinations with same momentum squared
    // The folders are created when running create_runs.sh. If they dont exist, 
    // a segmentation fault will occur
    // TODO: implement check for existence of folders

    for(int dirac = 0; dirac < number_of_dirac; ++dirac){
      for(int p = 0; p <= max_mom_squared; p++){
        for(int displ = 0; displ < number_of_displ; ++displ){

          sprintf(outfile, 
              "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d_unsuppressed/"
              "C2_pi+-_conf%04d.dat", 
              outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
              displ_min, displ_max, config_i);
          if((fp = fopen(outfile, "wb")) == NULL)
            std::cout << "fail to open outputfile: " << outfile << std::endl;

          for(int p_u = p_min; p_u < p_max; ++p_u){
            if(mom_squared[p_u] == p){

          		fwrite((double*) &(C2_mes[p_u][p_u][dirac][dirac]
                  [displ][displ][0]), 
                  sizeof(double), 2 * Lt, fp);
            }
          }

          fclose(fp);

        }
      }
    }

    for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
      for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
        for(int p1 = 0; p1 <= max_mom_squared; p1++){
         for(int p2 = p1; p2 <= max_mom_squared; p2++){
           for(int displ_u = 0; displ_u < number_of_displ; ++displ_u){
              for(int displ_d = 0; displ_d < number_of_displ; ++displ_d){

                printf("Writing to file: ");
                sprintf(outfile, 
                    "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
                    "C2_pi+-_conf%04d.dat", 
                    outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
                    p1, p2, displ_min + displ_u, displ_min + displ_d, config_i);
                printf("%s\n", outfile);
                if((fp = fopen(outfile, "wb")) == NULL)
                  std::cout << "fail to open outputfile" << outfile << std::endl;

                for(int p_u = p_min; p_u < p_max; ++p_u){
                  if(mom_squared[p_u] == p1){
                    for(int p_d = p_min; p_d < p_max; ++p_d){
                      if(mom_squared[p_d] == p2){

                        fwrite((double*) &(C2_mes[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d][0]), 
                            sizeof(double), 2 * Lt, fp);
                      }
                    }
                  }
                }

                fclose(fp);
  
              }
            }
          }
        }
      }
    }



    time = clock() - time;
    std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC << " seconds" << std::endl;

    //exit(0);

    // *************************************************************************
    // FOUR PT CONTRACTION 1 ***************************************************
    // *************************************************************************

    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
    time = clock();

    // displacement not supported for 4pt functions atm
    displ_min = 0;
    displ_max = 0;

    // setting the correlation function to zero

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u)
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = 
                  std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){     
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
             for(int p_u = p_min; p_u < p_max; ++p_u) {
               for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){
          
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
                                  [dirac_u][dirac_d][0][0]
                                  [t_source_1][t_sink_1][rnd1][rnd3]) *
                                (Corr[number_of_momenta - p_u - 1]
                                  [p_d][dirac_u][dirac_d][0][0]
                                  [t_source][t_sink][rnd2][rnd4]);
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

    // Normalization of 4pt-function. Accounts for all rnd-number combinations

    for(int t = 0; t < Lt; ++t){
      for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){     
        for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
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
    // see output to binary file for C2. 
    // write into folders with suffix "_unsuppressed". These only include
    // correlators of the diagonal matrix elements of the GEVP for which
    // the three-momentum remains unchanged for both quarks. Because the
    // quarks have to be back-to-back, for the offdiagonal elements this
    // cannot occur. The suppression can be interpreted as Zweig-suppressed
    // gluon exchange


    for(int dirac = 0; dirac < number_of_dirac; ++dirac){
      for(int p = 0; p <= max_mom_squared; p++){

        sprintf(outfile, 
            "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d_unsuppressed/"
            "C4_1_conf%04d.dat", 
            outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
            displ_min, displ_max, config_i);
        if((fp = fopen(outfile, "wb")) == NULL)
          std::cout << "fail to open outputfile" << std::endl;

        for(int p_u = p_min; p_u < p_max; ++p_u){
          if(mom_squared[p_u] == p){

            fwrite((double*) &(C4_mes[p_u][p_u][dirac][dirac][0]), 
                sizeof(double), 2 * Lt, fp);
          }
        }

        fclose(fp);

      }
    }

    // to build a GEVP, the correlators are written into a seperate folder
    // for every dirac structure, momentum, (entry of the GEVP matrix).
    // displacement is not supported at the moment

    for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
      for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
        for(int p1 = 0; p1 <= max_mom_squared; p1++){
          for(int p2 = p1; p2 <= max_mom_squared; p2++){

            sprintf(outfile, 
               "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
               "C4_1_conf%04d.dat", 
               outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
               p1, p2, displ_min, displ_max, config_i);
           if((fp = fopen(outfile, "wb")) == NULL)
             std::cout << "fail to open outputfile" << std::endl;

           for(int p_u = p_min; p_u < p_max; ++p_u){
              if(mom_squared[p_u] == p1){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(mom_squared[p_d] == p2){

                    fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), 
                        sizeof(double), 2 * Lt, fp);
                  }
                }
              }
            }

            fclose(fp);

          }
        }
      }
    }

    // output to terminal
//		printf("\n");
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
//		  printf("\tdirac    = %02d\n", dirac);
//      for(int offset = 0; offset <= max_mom_squared; offset++){
//        for(int p = 0; p <= max_mom_squared; p++){
//          if((p + offset) <= max_mom_squared){
//            for(int p_u = p_min; p_u < p_max; ++p_u){
//              if((rewr.mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
//                for(int p_d = p_min; p_d < p_max; ++p_d){
//                  if(rewr.mom_squared[p_d] == (p + offset)){
//            			  //printf(
//            				//  	"\t t\tRe(C4_1_con)\tIm(C4_1_con)\n\t----------------------------------\n");
////            			  for(int t1 = 0; t1 < Lt; ++t1){
////            				  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p_u][p_d][dirac][dirac][t1]),
////            				      imag(C4_mes[p_u][p_d][dirac][dirac][t1]));
////            			  }
//            			printf("\n");
//                  printf("p_u = %02d\tp_d = %02d\n", p_u, p_d);
//            		  }
//                }
//              }
//            }
//          }
//        printf("\n");
//        }
//      }
//    printf("\n");
//    }

    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

    // *************************************************************************
    // FOUR PT CONTRACTION 2 ***************************************************
    // *************************************************************************

    // setting the correlation function to zero
    std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
    time = clock();

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u)
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = 
                  std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt - 1; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){     
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
             for(int p_u = p_min; p_u < p_max; ++p_u) {
               for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(mom_squared[p_u] <= mom_squared[p_d]){

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
                                  [dirac_u][dirac_d][0][0][t_source_1][t_sink]
                                  [rnd1][rnd3]) *
                                (Corr[number_of_momenta - p_u - 1][p_d]
                                  [dirac_u][dirac_d][0][0][t_source][t_sink_1]
                                  [rnd2][rnd4]);
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

    // Normalization of 4pt-function. Accounts for all rnd-number combinations
    for(int t = 0; t < Lt; ++t){
      for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){     
        for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
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
    // see output to binary file for C2. 
    // write into folders with suffix "_unsuppressed". These only include
    // correlators of the diagonal matrix elements of the GEVP for which
    // the three-momentum remains unchanged for both quarks. Because the
    // quarks have to be back-to-back, for the offdiagonal elements this
    // cannot occur. The suppression can be interpreted as Zweig-suppressed
    // gluon exchange

    for(int dirac = 0; dirac < number_of_dirac; ++dirac){
      for(int p = 0; p <= max_mom_squared; p++){

        sprintf(outfile, 
            "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d_unsuppressed/"
            "C4_2_conf%04d.dat", 
            outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
            displ_min, displ_max, config_i);
        if((fp = fopen(outfile, "wb")) == NULL)
          std::cout << "fail to open outputfile" << std::endl;

        for(int p_u = p_min; p_u < p_max; ++p_u){
          if(mom_squared[p_u] == p){

            fwrite((double*) &(C4_mes[p_u][p_u][dirac][dirac][0]), 
                sizeof(double), 2 * Lt, fp);
          }
        }

        fclose(fp);

      }
    }

    // to build a GEVP, the correlators are written into a seperate folder
    // for every dirac structure, momentum, (entry of the GEVP matrix).
    // displacement is not supported at the moment

    for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
      for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
        for(int p1 = 0; p1 <= max_mom_squared; p1++){
          for(int p2 = p1; p2 <= max_mom_squared; p2++){

            sprintf(outfile, 
               "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
               "C4_2_conf%04d.dat", 
               outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
               p1, p2, displ_min, displ_max, config_i);
           if((fp = fopen(outfile, "wb")) == NULL)
             std::cout << "fail to open outputfile" << std::endl;

           for(int p_u = p_min; p_u < p_max; ++p_u){
              if(mom_squared[p_u] == p1){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(mom_squared[p_d] == p2){

                    fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), 
                        sizeof(double), 2 * Lt, fp);
                  }
                }
              }
            }

            fclose(fp);

          }
        }
      }
    }

    // output to terminal
//		printf("\n");
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
//		  printf("\tdirac    = %02d\n", dirac);
//      for(int offset = 0; offset <= max_mom_squared; offset++){
//        for(int p = 0; p <= max_mom_squared; p++){
//          if((p + offset) <= max_mom_squared){
//            for(int p_u = p_min; p_u < p_max; ++p_u){
//              if((rewr.mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
//                for(int p_d = p_min; p_d < p_max; ++p_d){
//                  if(rewr.mom_squared[p_d] == (p + offset)){
//            			  //printf(
//            				//  	"\t t\tRe(C4_2_con)\tIm(C4_2_con)\n\t----------------------------------\n");
////            			  for(int t1 = 0; t1 < Lt; ++t1){
////            				  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
////            				      imag(C4_mes[p][p][dirac][dirac][t1]));
////                    }
//            			  printf("\n");
//                    printf("p_u = %02d\tp_d = %02d\n", p_u, p_d);
//            		  }
//                }
//              }
//            }
//          }
//          printf("\n");
//        }
//      }
//      printf("\n");
//    }

    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

    // *************************************************************************
    // FOUR PT CONTRACTION 3 ***************************************************
    // *************************************************************************

    // TODO: check dirac indices. maybe dirac(t_source) and dirac(t_sink) have
    // to be equal or there may be four different structures rather than u- and
    // d-quark always having the same dirac structure
    // doesn't matter as long as all used dirac structures are equal

    std::cout << "\n\tcomputing the connected contribution of C4_3:\n";
    time = clock();

    // setting the correlation function to zero

    for(int p_u = 0; p_u < number_of_momenta; ++p_u)
      for(int p_d = 0; p_d < number_of_momenta; ++p_d)
        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u)
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d)
            for(int t1 = 0; t1 < Lt; ++t1)
              C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = 
                  std::complex<double>(0.0, 0.0);

    for(int t_source = 0; t_source < Lt; ++t_source){
      for(int t_sink = 0; t_sink < Lt; ++t_sink){

        int t_source_1 = (t_source + 1) % Lt;
        int t_sink_1 = (t_sink + 1) % Lt;

        // initialize basic->contraction[]
        // p_u = number_of_momenta/2 and the break; statement arrange
        // for one-to-all calculation in momentum space. (only one source
        // momentum is used. the first five are {(0,0,0), (0,0,1), 
        // (0,1,-1), (1,-1,-1), (0,0,2)}

        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
          for(int p = 0; p <= max_mom_squared; p++){
            for(int p_u = number_of_momenta/2; p_u < p_max; ++p_u){
              if(mom_squared[p_u] == p){
                basic->init_operator_u(0, t_source, t_sink, rewr, vdaggerv, 'b', p_u, 0);
                basic->init_operator_u(1, t_source_1, t_sink_1, rewr, vdaggerv, 'b', 
                                       number_of_momenta - p_u - 1, 0);
                break;
              }
            }
          }
        }

        // initialize basic->contraction_dagger[]
        // build all momenta for sinks

        for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
          for(int p = 0; p <= max_mom_squared; p++){
            for(int p_d = number_of_momenta/2; p_d > p_min; --p_d){
              if(mom_squared[p_d] == p){

                basic->init_operator_d(0, t_source_1, t_sink, rewr, vdaggerv, 'b', p_d, 0);
                basic->init_operator_d(1, t_source, t_sink_1, rewr, vdaggerv, 'b', 
                    number_of_momenta - p_d - 1, 0);
                break;
              }
            }
          }
        }

        // build 4pt-function C4_mes for pi^+pi^+ Equivalent two just summing
        // up the four-trace with same time difference between source and sink 
        // (all to all) for every dirac structure, momentum
        // displacement not supported at the moment
        // to build the trace with four matrices, build combinations 
        // X = D_d^-1(t_sink | t_source + 1) 
        //     Gamma D_u^-1(t_source + 1 | t_sink + 1) Gamma
        // Y = D_d^-1(t_sink + 1| t_source) 
        //     Gamma D_u^-1(t_source| t_sink) Gamma
        // these have dimension // (4 * quarks[0].number_of_dilution_E) x (4 * 
        //     quarks[0].number_of_dilution_E)
        // thus the multiplication in this order is fastest

        for(int dirac_1 = 0; dirac_1 < number_of_dirac; ++dirac_1){     
          for(int p = 0; p <= max_mom_squared; p++){
             for(int p_u = number_of_momenta / 2; p_u < p_max; ++p_u) {
              if(mom_squared[p_u] == p){
                for(int dirac_2 = 0; dirac_2 < number_of_dirac; ++dirac_2){
                  for(int p_d = p_min; p_d < p_max; ++p_d) {
                    if(p_u == p_d){

                      // initialisation of X. rnd loops and if-statements rule
                      // forbidden randomvector combinations (to improve 
                      // statistical error never use the same randomvector
                      // for different indices
    
                      basic->get_operator_g5(op_2, 0, dirac_ind.at(dirac_1), 
                          number_of_momenta - p_d - 1);
                      basic->get_operator_charged(op_3, 1, t_sink_1, rewr, 
                          dirac_ind.at(dirac_2), number_of_momenta - p_u - 1);
          
                      // second u quark: t_source_1 -> t_sink_1

                      for(int rnd3 = 1; rnd3 < number_of_rnd_vec; ++rnd3){
                        for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
                          if(rnd2 != rnd3){

                            // first d quark: t_sink_1 -> t_source

                            for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                              if(rnd4 != rnd3){

                                X[rnd3][rnd2][rnd4] = op_2[rnd3] * 
                                    op_3[rnd2][rnd4] ;
                              }
                            }
                          }
                        }
                      }

                      // initialisation of Y. see initialisation of X
    
                      basic->get_operator_g5(op_4, 1, dirac_ind.at(dirac_2), p_d);
                      basic->get_operator_charged(op_1, 0, t_sink, rewr, 
                          dirac_ind.at(dirac_1), p_u);
    
                      // first u quark: t_source -> t_sink

                      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                        for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){      

                          // second d quark: t_sink -> t_source_1

                          for(int rnd4 = 1; rnd4 < number_of_rnd_vec; ++rnd4){
                            if((rnd4 != rnd1) && (rnd4 != rnd3)){

                              Y[rnd4][rnd1][rnd3] = op_4[rnd4] * 
                                  op_1[rnd1][rnd3];
                            }
                          }
                        }
                      }
              
                      // complete diagramm. combine X and Y to four-trace
                      // C4_mes = tr(D_u^-1(t_source| t_sink) Gamma 
                      //     D_d^-1(t_sink | t_source + 1) Gamma 
                      //     D_u^-1(t_source + 1 | t_sink + 1) Gamma 
                      //     D_d^-1(t_sink + 1| t_source) Gamma)
                      // every quark line must have its own random vec

                      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                        for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
                          for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){      
                            if((rnd2 != rnd1) && (rnd2 != rnd3)){
                              for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
                                if((rnd4 != rnd1) && (rnd4 != rnd3)){

                                  C4_mes[p_u][p_d][dirac_1][dirac_2]
                                      [abs((t_sink - t_source - Lt) % Lt)] += 
                                    ((X[rnd3][rnd2][rnd4] * 
                                      Y[rnd4][rnd1][rnd3]).trace());
                                }
                              }
                            }
                          }
                        }
                      }
    
                    }
                  }
                }

                break;

              }
            }
          }
        }

      }
    }

    // Normalization of 4pt-function. Accounts for all rnd-number combinations

    for(int p1 = 0; p1 < number_of_momenta; ++p1)
      for(int p2 = 0; p2 < number_of_momenta; ++p2)
        for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u)
          for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d)
            for(int t = 0; t < Lt; ++t)
              C4_mes[p1][p2][dirac_u][dirac_d][t] /= norm1;


    // output to binary file

    // see output to binary file for C2. 
    // for the C4_3 diagram the four propagators are connected in the same
    // trace. Thus there are no gluon lines which could be cut to create a
    // disconnected diagrams and thus no Zweig suppression.
    // To build a GEVP, the correlators are written into a seperate folder
    // for every dirac structure, momentum, (entry of the GEVP matrix).
    // displacement is not supported at the moment

    for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
      for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
        for(int p1 = 0; p1 <= max_mom_squared; p1++){
          for(int p2 = p1; p2 <= max_mom_squared; p2++){

            sprintf(outfile, 
                "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
                "C4_3_conf%04d.dat", 
                outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
                p1, p2, displ_min, displ_max, config_i);
            if((fp = fopen(outfile, "wb")) == NULL)
              std::cout << "fail to open outputfile" << std::endl;

            for(int p_u = number_of_momenta / 2; p_u < p_max; ++p_u){
              if(mom_squared[p_u] == p1){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(mom_squared[p_d] == p2){

                    fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), 
                        sizeof(double), 2 * Lt, fp);
                  }
                }

                break;

              }
            }

            fclose(fp);

          }
        }
      }
    }

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
//    printf("\n");
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
//      printf("\tdirac    = %02d\n", dirac);
//      for(int p = p_min; p < p_max; ++p) {
//        printf("\tmomentum = %02d\n", p);
//        //printf(
//        //    "\t t\tRe(C4_3_con)\tIm(C4_3_con)\n\t----------------------------------\n");
//        for(int t1 = 0; t1 < Lt; ++t1){
//          printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
//              imag(C4_mes[p][p][dirac][dirac][t1]));
//        }
//        printf("\n");
//      }
//      printf("\n");
//    }


    time = clock() - time;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

    // *************************************************************************
    // FOUR PT CONTRACTION 4 ***************************************************
    // *************************************************************************

    // identical to FOUR PT CONTRACTION 3

  } // loop over configs ends here

  // TODO: freeing all memory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  delete basic;
}

