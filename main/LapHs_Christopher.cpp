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

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "boost/multi_array.hpp"

#include "GlobalData.h"
#include "BasicOperator.h"
#include "ReadWrite.h"

typedef boost::multi_array<Eigen::MatrixXcd, 2> dim2_eigen_array;
typedef boost::multi_array<Eigen::MatrixXcd, 3> dim3_eigen_array;
typedef std::complex<double> cmplx;
typedef boost::multi_array<cmplx, 4> dim4_array;
typedef boost::multi_array<cmplx, 5> dim5_array;
typedef boost::multi_array<cmplx, 7> dim7_array;
typedef boost::multi_array<cmplx, 10> dim10_array;

int main (int ac, char* av[]) {

	//Eigen::initParallel();

	// Reading in global parameters from input file
	GlobalData* global_data = GlobalData::Instance();
	global_data->read_parameters(ac, av);

	//Eigen::setNbThreads(4);

	// global variables from input file needed in main function
	const int Lt = global_data->get_Lt();
	const int end_config = global_data->get_end_config();
	const int delta_config = global_data->get_delta_config();
	const int start_config = global_data->get_start_config();
	const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_squared = number_of_max_mom * number_of_max_mom;

	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
	//const int number_of_inversions = quarks[0].number_of_dilution_T
  //			* quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;

	clock_t time;

	const std::complex<double> I(0.0, 1.0);

	char outfile[400];
	FILE *fp = NULL;

  // everything for Input files
  ReadWrite rewr = ReadWrite();
//  ReadWrite* rewr = &r;

	// everything for operator handling
	BasicOperator basic = BasicOperator(&rewr);
//  BasicOperator* basic = &b;

  // ***************************************************************************
	// ***************************************************************************
	// dirac indices and momenta to calculate ************************************ 
	// ***************************************************************************
	// ***************************************************************************

  //TODO: put that into infile and get it in build_source matrix from globaldata
  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();
  const int number_of_dirac = dirac_max - dirac_min + 1;

  const int displ_min = global_data->get_displ_min();
  const int displ_max = global_data->get_displ_max();
  const int number_of_displ = displ_max - displ_min + 1;

  const int p_min = 0; //rewr.number_of_momenta/2;
  const int p_max = rewr.number_of_momenta;

  std::string outpath = global_data->get_output_path() + "/" + global_data->
      get_name_lattice();

	// ***************************************************************************
	// memory allocation *********************************************************
	// memory for the operator in Dirac and eigenvector space and in time
	// ***************************************************************************
  // ***************************************************************************

  // abbreviations for clearer memory allocation. Wont be used in loops and 
  //when building the contractions

  const size_t nmom = rewr.number_of_momenta;
  const size_t nrnd = number_of_rnd_vec;
  const size_t ndir = number_of_dirac;
  const size_t ndis = number_of_displ;

	// memory for the correlation function

  dim7_array C2_mes(boost::extents[nmom][nmom][ndir][ndir][ndis][ndis][Lt]);
  //TODO: dont need the memory for p_u^2 > p_d^2
  dim10_array Corr(boost::extents[nmom][nmom][ndir][ndir][ndis][ndis][Lt][Lt][nrnd][nrnd]);

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
  // 1 -> u-quark; 2 -> d-quarks
  dim2_eigen_array op_1(boost::extents[nrnd][nrnd]);
  std::vector<Eigen::MatrixXcd> op_2(number_of_rnd_vec);

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

    // initialize everything regarding I/O. 
    // read_lime_gauge_field_doubleprec_timeslices() can be skipped
    // when no displacemnent is used (displ_min = displ_max = 0)

		rewr.read_perambulators_from_file(config_i);
		rewr.read_rnd_vectors_from_file(config_i);
    rewr.read_lime_gauge_field_doubleprec_timeslices(config_i);
    rewr.build_source_matrix(config_i, p_min, p_max);


		// *************************************************************************
		// TWO PT CONTRACTION 1 ****************************************************
		// *************************************************************************


		std::cout << "\n\tcomputing the connected contribution of pi_+/-:\n";
		time = clock();

		// setting the correlation function to zero

		for(int p1 = 0; p1 < rewr.number_of_momenta; ++p1)
		  for(int p2 = 0; p2 < rewr.number_of_momenta; ++p2)
			  for(int dirac1 = 0; dirac1 < number_of_dirac; ++dirac1)
			    for(int dirac2 = 0; dirac2 < number_of_dirac; ++dirac2)
				    for(int displ1 = 0; displ1 < number_of_displ; ++displ1)
			        for(int displ2 = 0; displ2 < number_of_displ; ++displ2)
      			    for(int t1 = 0; t1 < Lt; ++t1)
           				for(int t1 = 0; t1 < Lt; ++t1)
           					C2_mes[p1][p2][dirac1][dirac2][displ1][displ2][t1] = std::complex<double>(0.0, 0.0);

		for(int p1 = 0; p1 < rewr.number_of_momenta; ++p1)
		  for(int p2 = 0; p2 < rewr.number_of_momenta; ++p2)
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

    // initializing of Corr: calculate all two-operator traces of the form tr(u \Gamma \bar{d})
    // build all combinations of momenta, dirac_structures and displacements as specified in
    // infile

    for(int displ_u = 0; displ_u < number_of_displ; displ_u++){
      for(int displ_d = 0; displ_d < number_of_displ; displ_d++){

    		  for(int t_source = 0; t_source < Lt; ++t_source){
    		  	for(int t_sink = 0; t_sink < Lt; ++t_sink){
    
    	        for(int p = p_min; p < p_max; ++p) {
                // initialize contraction[rnd_i] = perambulator * basicoperator
                // = D_u^-1
                // choose 'i' for interlace or 'b' for block time dilution scheme
                // TODO: get that from input file
                // choose 'c' for charged or 'u' for uncharged particles
                basic.init_operator_u(0, t_source, t_sink, &rewr, 'b', p, displ_min + displ_u);
                basic.init_operator_d(0, t_source, t_sink, &rewr, 'b', p, displ_min + displ_d);
              }
    
              for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
    	          for(int p_u = p_min; p_u < p_max; ++p_u) {
                  // code for pi+-
    
    	            // "multiply contraction[rnd_i] with gamma structure"
                  // contraction[rnd_i] are the columns of D_u^-1 which get
                  // reordered by gamma multiplication. No actual multiplication
                  // is carried out
                  basic.get_operator_charged(op_1, 0, t_sink, &rewr, dirac_min + dirac_u, p_u);
    
                  for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
                    for(int p_d = p_min; p_d < p_max; ++p_d) {
                      if(rewr.mom_squared[p_u] <= rewr.mom_squared[p_d]){
      
                        // same as get_operator but with gamma_5 trick. D_u^-1 is
                        // daggered and multipied with gamma_5 from left and right
                        // the momentum is changed to reflect the switched sign in
                        // the momentum exponential for pi_+-
                        basic.get_operator_g5(op_2, 0, dirac_min + dirac_d, p_d);
           
                        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                          for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){

                            // build all 2pt traces leading to C2_mes
                            // Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)
                            
                            Corr[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d][t_source][t_sink][rnd1][rnd2] = 
                              (op_2[rnd2] * op_1[rnd1][rnd2]).trace();
          
          
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

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

   	    for(int p_u = p_min; p_u < p_max; ++p_u) {
   	      for(int p_d = p_min; p_d < p_max; ++p_d) {
            if(rewr.mom_squared[p_u] <= rewr.mom_squared[p_d]){
              for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
                for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
                  for(int displ_u = 0; displ_u < number_of_displ; displ_u++){
                    for(int displ_d = 0; displ_d < number_of_displ; displ_d++){
            
                      for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                        for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
              
                          // building Correlation function C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                          // TODO: find signflip of imaginary part
                          // TODO: is C2_mes[dirac][p] better?

                          C2_mes[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d][abs((t_sink - t_source - Lt) % Lt)] += 
                              Corr[p_u][rewr.number_of_momenta - p_d - 1][dirac_u][dirac_d]
                              [displ_u][displ_d][t_source][t_sink][rnd1][rnd2];
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

    // normalization of 2pt-function. Accounts for all rnd number combinations

		double norm3 = Lt * number_of_rnd_vec * (number_of_rnd_vec - 1) * 0.5;
    for(int p_u = p_min; p_u < p_max; ++p_u) {
      for(int p_d = p_min; p_d < p_max; ++p_d) {
        if(rewr.mom_squared[p_u] <= rewr.mom_squared[p_d]){
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


    // output to binary file
    // to build a GEVP, the correlators are written into a seperate folder
    // for every dirac structure, momentum, displacement (entry of the GEVP
    // matrix. In the folders a file is created for every configuration which
    // contains all momentum combinations with same momentum squared
    // @Christopher: you can adjust these as you need.
    // The folders are created when running create_runs.sh. If they dont exist, 
    // a segmentation fault will occur
    // TODO: implement check for existence of folders

    for(int dirac_u = 0; dirac_u < number_of_dirac; ++dirac_u){
      for(int dirac_d = 0; dirac_d < number_of_dirac; ++dirac_d){
        for(int p1 = 0; p1 <= max_mom_squared; p1++){
         for(int p2 = p1; p2 <= max_mom_squared; p2++){
           for(int displ_u = 0; displ_u < number_of_displ; ++displ_u){
              for(int displ_d = 0; displ_d < number_of_displ; ++displ_d){

		            sprintf(outfile, 
                    "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/C2_pi+-_conf%04d.dat", 
                    outpath.c_str(), dirac_min + dirac_u, dirac_max + dirac_d, 
                    p1, p2, displ_u, displ_d, config_i);
		            if((fp = fopen(outfile, "wb")) == NULL)
		            	std::cout << "fail to open outputfile" << std::endl;

                for(int p_u = p_min; p_u < p_max; ++p_u){
                  if(rewr.mom_squared[p_u] == p1){
                    for(int p_d = p_min; p_d < p_max; ++p_d){
                      if(rewr.mom_squared[p_d] == p2){

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

    // output to terminal
//		printf("\n");
//    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
//		  printf("\tdirac    = %02d\n", dirac);
//      for(int p = 0; p <= max_mom_squared; p++){
//        printf("\tmomentum_u = %02d\n", p);
//        printf("\tmomentum_d = %02d\n", p);
//        for(int p_u = p_min; p_u < p_max; ++p_u){
//          if((rewr.mom_squared[p_u] == p)){
//      			//printf(
//      			//  	"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
////            		  for(int t1 = 0; t1 < Lt; ++t1){
////            			  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
////            			      imag(C2_mes[p_u][p_u][dirac][t1]));
////            		  }
//            printf("\n");
//            printf("p_u = %02d\n", p_u);
//          }
//        }
//      }
//
//      for(int p = 1; p <= max_mom_squared; p++){
//        printf("\tmomentum_u = %02d\n", 0);
//        printf("\tmomentum_d = %02d\n", p);
//        for(int p_u = p_min; p_u < p_max; ++p_u){
//          if((rewr.mom_squared[p_u] == p)){
//      			//printf(
//      			//  	"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
////            		  for(int t1 = 0; t1 < Lt; ++t1){
////            			  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
////            			      imag(C2_mes[p_u][p_u][dirac][t1]));
////            		  }
//            printf("\n");
//            printf("p_u = %02d\n", p_u);
//          }
//        }
//      }
//
//    }

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

	} // loop over configs ends here

	return 0;
}

