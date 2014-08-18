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
typedef boost::multi_array<cmplx, 8> dim8_array;

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
	BasicOperator basic = BasicOperator(rewr);
//  BasicOperator* basic = &b;

  // ***************************************************************************
	// ***************************************************************************
	// dirac indices and momenta to calculate ************************************ 
	// ***************************************************************************
	// ***************************************************************************

  //TODO: put that into infile and get it in build_source matrix from globaldata
  const int dirac_min = 0; //global_data->get_dirac_min();
  const int dirac_max = 0; //global_data->get_dirac_max();

  const int displ_min = 0; //global_data->get_displ_min()?;
  const int displ_max = 0;

  const int p_min = 0; //rewr->number_of_momenta/2;
  const int p_max = rewr->number_of_momenta;

  std::string outpath = global_data->get_output_path() + "/" + global_data->
      get_name_lattice();

	// ***************************************************************************
	// memory allocation *********************************************************
	// memory for the operator in Dirac and eigenvector space and in time
	// ***************************************************************************
  // ***************************************************************************

	Eigen::MatrixXcd op_D_tsource = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsource_1 = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsink = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,

			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsink_1 = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);

  dim3_eigen_array X(boost::extents[number_of_rnd_vec][number_of_rnd_vec][number_of_rnd_vec]);
  dim3_eigen_array Y(boost::extents[number_of_rnd_vec][number_of_rnd_vec][number_of_rnd_vec]);
	for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
		for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++){
      for(int rnd3 = 0; rnd3 < number_of_rnd_vec; rnd3++){
  			X[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
  					4 * quarks[0].number_of_dilution_E, 4 * quarks[0].number_of_dilution_E);
  			Y[rnd1][rnd2][rnd3] = Eigen::MatrixXcd::Zero(
  					4 * quarks[0].number_of_dilution_E, 4 * quarks[0].number_of_dilution_E);
      }
		}
	}

	// memory for the correlation function
  const size_t nmom = rewr->number_of_momenta;
  const size_t nrnd = number_of_rnd_vec;
  dim4_array C2_mes(boost::extents[nmom][nmom][1][Lt]);
  dim4_array C2_dis(boost::extents[nmom][nmom][1][Lt]);
  dim5_array C4_mes(boost::extents[nmom][nmom][1][1][Lt]);
  dim8_array Corr(boost::extents[nmom][nmom][1][1][Lt][Lt][nrnd][nrnd]);

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
  dim2_eigen_array op_1(boost::extents[nrnd][nrnd]);
  dim2_eigen_array op_3(boost::extents[nrnd][nrnd]);
  std::vector<Eigen::MatrixXcd> op_2(number_of_rnd_vec);
  std::vector<Eigen::MatrixXcd> op_4(number_of_rnd_vec);
  std::vector<Eigen::MatrixXcd> op_5(number_of_rnd_vec);
  std::vector<Eigen::MatrixXcd> op_6(number_of_rnd_vec);

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

    op_5[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
        4 * number_of_eigen_vec);
    op_6[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
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

		rewr->read_perambulators_from_file(config_i);
//		rewr->read_eigenvectors_from_file(config_i);
		rewr->read_rnd_vectors_from_file(config_i);
//    rewr->read_lime_gauge_field_doubleprec_timeslices(config_i);
    rewr->build_source_matrix(config_i, p_min, p_max);


		// *************************************************************************
		// TWO PT CONTRACTION 1 ****************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of pi_+/-:\n";
		time = clock();

		for(int p_u = 0; p_u < rewr->number_of_momenta; ++p_u)
		  for(int p_d = 0; p_d < rewr->number_of_momenta; ++p_d)
  			for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
  				for(int t1 = 0; t1 < Lt; ++t1)
  					C2_mes[p_u][p_d][dirac][t1] = std::complex<double>(0.0, 0.0);

		for(int p1 = 0; p1 < rewr->number_of_momenta; ++p1)
		  for(int p2 = 0; p2 < rewr->number_of_momenta; ++p2)
			  for(int dirac1 = dirac_min; dirac1 < dirac_max + 1; ++dirac1)
			    for(int dirac2 = dirac_min; dirac2 < dirac_max + 1; ++dirac2)
				    for(int t1 = 0; t1 < Lt; ++t1)
				      for(int t2 = 0; t2 < Lt; ++t2)
                for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++)
                  for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++)
                    Corr[p1][p2][dirac1][dirac2][t1][t2][rnd1][rnd2] = 
                        std::complex<double>(0.0, 0.0);


#if 0
		for(int p = 0; p < rewr->number_of_momenta; ++p)
			for(int dirac = 0; dirac < 16; ++dirac)
        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++)
  				for(int t1 = 0; t1 < Lt; ++t1)
  					C2_dis[p][dirac][rnd1][t1] = std::complex<double>(0.0, 0.0);

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

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

	      for(int p = p_min; p < p_max; ++p) {
          // initialize contraction[rnd_i] = perambulator * basicoperator
          // = D_u^-1
          // choose 'i' for interlace or 'b' for block time dilution scheme
          // choose 'c' for charged or 'u' for uncharged particles
          basic->init_operator_u(0, t_source, t_sink, rewr, 'b', p, 0);
          basic->init_operator_d(0, t_source, t_sink, rewr, 'b', p, 0);
        }

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
	        for(int p_u = p_min; p_u < p_max; ++p_u) {
            // code for pi+-

	          // "multiply contraction[rnd_i] with gamma structure"
            // contraction[rnd_i] are the columns of D_u^-1 which get
            // reordered by gamma multiplication. No actual multiplication
            // is carried out
            basic->get_operator_charged(op_1, 0, t_sink, rewr, 5, p_u);

            for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
              for(int p_d = p_min; p_d < p_max; ++p_d) {
                if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){
  
                  // same as get_operator but with gamma_5 trick. D_u^-1 is
                  // daggered and multipied with gamma_5 from left and right
                  // the momentum is changed to reflect the switched sign in
                  // the momentum exponential for pi_+-
                  basic->get_operator_g5(op_2, 0, 5, p_d);
     
                  for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                    for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
                        // building Correlation function get quite intuitive
                        // C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                        // TODO: find signflip of imaginary part
                        // TODO: is C2_mes[dirac][p] better?
                      Corr[p_u][p_d][dirac_u][dirac_d][t_source][t_sink][rnd1][rnd2] = 
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



		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

        for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
    	    for(int p_u = p_min; p_u < p_max; ++p_u) {
    	      for(int p_d = p_min; p_d < p_max; ++p_d) {
    
              for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
                for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
      
                  C2_mes[p_u][p_d][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
                      Corr[p_u][rewr->number_of_momenta - p_d - 1][dirac][dirac]
                      [t_source][t_sink][rnd1][rnd2];
                }
              }
      
            }
          }
        }

      }
    }

#if 0
// code for pi0 

            basic->init_operator_uncharged(t_source, t_sink, rewr, 'b', p, 0);
            basic->get_operator_uncharged(op_5, dirac);
    
            basic->init_operator_uncharged(t_sink, t_source, rewr, 'b', 
                rewr->number_of_momenta - p - 1, 0);
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
		sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_pi+-_conf%04d.dat", 
        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, 
        displ_min, displ_max, config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
      for(int p = 0; p <= max_mom_squared; p++){
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if(rewr->mom_squared[p_u] == p){
		        fwrite((double*) &(C2_mes[p_u][p_u][dirac][0]), sizeof(double), 2 * Lt, fp);
          }
        }
      }

      for(int p = 1; p <= max_mom_squared; p++){
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if(rewr->mom_squared[p_u] == p){
			      fwrite((double*) &(C2_mes[rewr->number_of_momenta / 2][p_u][dirac][0]), sizeof(double), 2 * Lt, fp);
          }
        }
      }

    }
		fclose(fp);

#if 0
		sprintf(outfile, 
        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_pi+-_conf%04d.dat", 
        outpath.c_str(), dirac_min, dirac_max, 0, 
        displ_min, displ_max, config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
	  for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
		  fwrite((double*) C2_mes[rewr->number_of_momenta/2]
          [rewr->number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
		fclose(fp);

    for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
  		sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C2_dis_u_rnd%02d_conf%04d.dat", 
          outpath.c_str(), dirac_min, dirac_max, number_of_max_mom, displ_min, 
          displ_max, rnd_i, config_i);
  		if((fp = fopen(outfile, "wb")) == NULL)
  			std::cout << "fail to open outputfile" << std::endl;
  		for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
  		  for(int p = 0; p < rewr->number_of_momenta; ++p)
  				fwrite((double*) C2_dis[p][dirac][rnd_i], sizeof(double), 2 * Lt, fp);
  		fclose(fp);
    }
#endif


    // output to terminal
		printf("\n");
    for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
		  printf("\tdirac    = %02d\n", dirac);
      for(int p = 0; p <= max_mom_squared; p++){
        printf("\tmomentum_u = %02d\n", p);
        printf("\tmomentum_d = %02d\n", p);
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if((rewr->mom_squared[p_u] == p)){
      			//printf(
      			//  	"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
//            		  for(int t1 = 0; t1 < Lt; ++t1){
//            			  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
//            			      imag(C2_mes[p_u][p_u][dirac][t1]));
//            		  }
            printf("\n");
            printf("p_u = %02d\n", p_u);
          }
        }
      }

      for(int p = 1; p <= max_mom_squared; p++){
        printf("\tmomentum_u = %02d\n", 0);
        printf("\tmomentum_d = %02d\n", p);
        for(int p_u = p_min; p_u < p_max; ++p_u){
          if((rewr->mom_squared[p_u] == p)){
      			//printf(
      			//  	"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
//            		  for(int t1 = 0; t1 < Lt; ++t1){
//            			  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[p_u][p_u][dirac][t1]),
//            			      imag(C2_mes[p_u][p_u][dirac][t1]));
//            		  }
            printf("\n");
            printf("p_u = %02d\n", p_u);
          }
        }
      }

    }

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);


		// *************************************************************************
		// FOUR PT CONTRACTION 1 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
		time = clock();

		for(int p_u = 0; p_u < rewr->number_of_momenta; ++p_u)
		  for(int p_d = 0; p_d < rewr->number_of_momenta; ++p_d)
		  	for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u)
		  	  for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d)
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
                if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){

  
      //            basic->init_operator_charged(t_source_1, t_sink_1, rewr, 'b', p, 0);
      //    
      //            basic->get_operator_charged(op_1, rewr, dirac, t_sink_1);
      //            basic->get_operator_g5(op_2, dirac);
      //    
      //            basic->init_operator_charged(t_source, t_sink, rewr, 'b', 
      //                rewr->number_of_momenta - p - 1, 0);
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
                                  [rewr->number_of_momenta - p_d - 1]
                                  [dirac_u][dirac_d]
                                  [t_source_1][t_sink_1][rnd1][rnd3]) *
                                (Corr[rewr->number_of_momenta - p_u - 1]
                                  [p_d][dirac_u][dirac_d]
                                  [t_source][t_sink][rnd2][rnd4]);
      //    							      C2_mes[p][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
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
              if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){
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
            if(rewr->mom_squared[p_u] == p){
          		fwrite((double*) &(C4_mes[p_u][p_u][dirac_u][dirac_d][0]), sizeof(double), 2 * Lt, fp);
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
                if(rewr->mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(rewr->mom_squared[p_d] == (p + offset)){
            			    fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), sizeof(double), 2 * Lt, fp);
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

//		sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_1_conf%04d.dat", 
//        outpath.c_str(), dirac_min, dirac_max, 0, displ_min, 
//        displ_max, config_i);
//		if((fp = fopen(outfile, "wb")) == NULL)
//			std::cout << "fail to open outputfile" << std::endl;
//		for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
//			fwrite((double*) C2_mes[rewr->number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
//		fclose(fp);

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
              if((rewr->mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(rewr->mom_squared[p_d] == (p + offset)){
            			  //printf(
            				//  	"\t t\tRe(C4_1_con)\tIm(C4_1_con)\n\t----------------------------------\n");
//            			  for(int t1 = 0; t1 < Lt; ++t1){
//            				  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p_u][p_d][dirac][dirac][t1]),
//            				      imag(C4_mes[p_u][p_d][dirac][dirac][t1]));
//            			  }
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

		// *************************************************************************
		// FOUR PT CONTRACTION 2 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
		time = clock();

		for(int p_u = 0; p_u < rewr->number_of_momenta; ++p_u)
		  for(int p_d = 0; p_d < rewr->number_of_momenta; ++p_d)
		  	for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u)
		  	  for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d)
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
                if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){


      //            basic->init_operator_charged(t_source_1, t_sink, rewr, 'b', p, 0);
      //    
      //            basic->get_operator_charged(op_1, rewr, dirac, t_sink);
      //            basic->get_operator_g5(op_2, dirac);
      //    
      //            basic->init_operator_charged(t_source, t_sink_1, rewr, 'b', 
      //                rewr->number_of_momenta - p - 1, 0);
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
                                (Corr[p_u][rewr->number_of_momenta - p_d - 1]
                                  [dirac_u][dirac_d][t_source_1][t_sink][rnd1][rnd3]) *
                                (Corr[rewr->number_of_momenta - p_u - 1][p_d]
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
              if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){
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
            if(rewr->mom_squared[p_u] == p){
              fwrite((double*) &(C4_mes[p_u][p_u][dirac_u][dirac_d][0]), sizeof(double), 2 * Lt, fp);
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
                if(rewr->mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(rewr->mom_squared[p_d] == (p + offset)){
                      fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), sizeof(double), 2 * Lt, fp);
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

//		sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_2_conf%04d.dat", 
//        outpath.c_str(), dirac_min, dirac_max, 0, displ_min, 
//        displ_max, config_i);
//		if((fp = fopen(outfile, "wb")) == NULL)
//			std::cout << "fail to open outputfile" << std::endl;
//		for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
//			fwrite((double*) C2_mes[rewr->number_of_momenta/2][dirac], sizeof(double), 2 * Lt, fp);
//		fclose(fp);

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
              if((rewr->mom_squared[p_u] == p) && ((p + offset) <= max_mom_squared)){
                for(int p_d = p_min; p_d < p_max; ++p_d){
                  if(rewr->mom_squared[p_d] == (p + offset)){
            			  //printf(
            				//  	"\t t\tRe(C4_2_con)\tIm(C4_2_con)\n\t----------------------------------\n");
//            			  for(int t1 = 0; t1 < Lt; ++t1){
//            				  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
//            				      imag(C4_mes[p][p][dirac][dirac][t1]));
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


		// *************************************************************************
		// FOUR PT CONTRACTION 3 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_3:\n";
		time = clock();
		for(int p_u = 0; p_u < rewr->number_of_momenta; ++p_u)
  		for(int p_d = 0; p_d < rewr->number_of_momenta; ++p_d)
  			for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u)
  			  for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d)
  			  	for(int t1 = 0; t1 < Lt; ++t1)
  			  		C4_mes[p_u][p_d][dirac_u][dirac_d][t1] = std::complex<double>(0.0, 0.0);

		for(int t_source = 0; t_source < Lt; ++t_source){
      std::cout << "\tt_source = " << t_source << std::endl;
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

				int t_source_1 = (t_source + 1) % Lt;
				int t_sink_1 = (t_sink + 1) % Lt;



        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
          for(int p = 0; p <= max_mom_squared; p++){
            for(int p_u = rewr->number_of_momenta/2; p_u < p_max; ++p_u){
              if(rewr->mom_squared[p_u] == p){
                basic->init_operator_u(0, t_source, t_sink, rewr, 'b', p_u, 0);
                basic->init_operator_u(1, t_source_1, t_sink_1, rewr, 'b', 
                    rewr->number_of_momenta - p_u - 1, 0);
//                std::cout << "p_u = " << p_u << std::endl;
                break;
              }
            }
          }
        }


        for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
          for(int p_d = p_min; p_d < p_max; ++p_d){
            basic->init_operator_d(0, t_source_1, t_sink, rewr, 'b', p_d, 0);
            basic->init_operator_d(1, t_source, t_sink_1, rewr, 'b', 
                rewr->number_of_momenta - p_d - 1, 0);
          }
        }

        for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){     
          for(int p = 0; p <= max_mom_squared; p++){
   	        for(int p_u = rewr->number_of_momenta / 2; p_u < p_max; ++p_u) {
              if(rewr->mom_squared[p_u] == p){
                for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
         	        for(int p_d = p_min; p_d < p_max; ++p_d) {
                    if(rewr->mom_squared[p_u] <= rewr->mom_squared[p_d]){
    
                      basic->get_operator_g5(op_2, 0, 5, 
                          rewr->number_of_momenta - p_d - 1);
                      basic->get_operator_charged(op_3, 1, t_sink_1, rewr, 5, 
                          rewr->number_of_momenta - p_u - 1);
          
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
    
                      basic->get_operator_g5(op_4, 1, 5, p_d);
                      basic->get_operator_charged(op_1, 0, t_sink, rewr, 5, p_u);
    
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
                                      ((X[rnd3][rnd2][rnd4] * Y[rnd4][rnd1][rnd3]).trace());
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

		for(int t = 0; t < Lt; ++t)
			for(int p = 0; p < rewr->number_of_momenta; ++p)
			  for(int p1 = 0; p1 < rewr->number_of_momenta; ++p1)
				for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac)
					C4_mes[p][p1][dirac][dirac][t] /= norm1;


    // output to binary file
//		sprintf(outfile, 
//        "%s/dirac_%02d_%02d_p_0_%01d_displ_%01d_%01d/C4_3_conf%04d_diag.dat", 
//        outpath.c_str(), dirac_min, dirac_max, max_mom_squared, displ_min, 
//        displ_max, config_i);
//		if((fp = fopen(outfile, "wb")) == NULL)
//			std::cout << "fail to open outputfile" << std::endl;
//    for(int dirac_u = dirac_min; dirac_u < dirac_max + 1; ++dirac_u){
//      for(int dirac_d = dirac_min; dirac_d < dirac_max + 1; ++dirac_d){
//        for(int p = 0; p <= max_mom_squared; p++){
//          for(int p_u = p_min; p_u < p_max; ++p_u){
//            if(rewr->mom_squared[p_u] == p){
//              fwrite((double*) C4_mes[p_u][p_u][dirac_u][dirac_d], sizeof(double), 2 * Lt, fp);
//            }
//          }
//        }
//      }
//    }
//		fclose(fp);

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
              for(int p_u = rewr->number_of_momenta / 2; p_u < p_max; ++p_u){
                if(rewr->mom_squared[p_u] == p){
                  for(int p_d = p_min; p_d < p_max; ++p_d){
                    if(rewr->mom_squared[p_d] == (p + offset)){
//                      std::cout << "p_u = " << p_u << "\tp_d = " << p_d << std::endl;
				              fwrite((double*) &(C4_mes[p_u][p_d][dirac_u][dirac_d][0]), sizeof(double), 2 * Lt, fp);
                    }
                  }
                  break;
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
			fwrite((double*) C4_mes[rewr->number_of_momenta/2]
          [rewr->number_of_momenta/2][dirac][dirac], sizeof(double), 2 * Lt, fp);
		fclose(fp);
#endif

    // output to terminal
//		printf("\n");
//		for(int dirac = dirac_min; dirac < dirac_max + 1; ++dirac){
//			printf("\tdirac    = %02d\n", dirac);
//      for(int p = p_min; p < p_max; ++p) {
//        printf("\tmomentum = %02d\n", p);
//			  //printf(
//				//  	"\t t\tRe(C4_3_con)\tIm(C4_3_con)\n\t----------------------------------\n");
//			  for(int t1 = 0; t1 < Lt; ++t1){
//				  printf("\t%02d\t%.5e\t%.5e\n", t1, real(C4_mes[p][p][dirac][dirac][t1]),
//				      imag(C4_mes[p][p][dirac][dirac][t1]));
//			  }
//			  printf("\n");
//		  }
//		  printf("\n");
//    }


    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
		
		// *************************************************************************
		// FOUR PT CONTRACTION 4 ***************************************************
		// *************************************************************************

    // identical to FOUR PT CONTRACTION 3

	} // loop over configs ends here

	// TODO: freeing all memory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//delete basic;
	return 0;
}

