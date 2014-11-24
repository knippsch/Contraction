#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::compute_meson_small_traces(const int t_source, 
                                                   const int t_sink){

  const size_t Lt = global_data->get_Lt();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const int dilT = quarks[0].number_of_dilution_T;
  // TODO: must be changed in GlobalData {
  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

  for(int displ_u = 0; displ_u < nb_dis; displ_u++){
  for(int displ_d = 0; displ_d < nb_dis; displ_d++){ 
    for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){
    for(int p_u = 0; p_u < nb_mom; ++p_u) {
      for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p_d = 0; p_d < nb_mom; ++p_d) {
        // TODO: A collpase of both random vectors might be better but
        //       must be done by hand because rnd2 starts from rnd1+1
        #pragma omp parallel for schedule(dynamic)
        for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
        for(int rnd2 = rnd1+1; rnd2 < nb_rnd; ++rnd2){
          // build all 2pt traces leading to C2_mes
          // Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)
          // TODO: Just a workaround
          std::array<double, 4> bla = {{1., 1., -1., -1.}};
          for(size_t block = 0; block < 4; block++){
            // TODO: dilution scheme in time should be choosable
            Corr[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d]
              [t_source][t_sink][rnd1][rnd2] += bla[block] *
              ((basic.get_operator(t_source, t_sink/dilT, 1, dirac_u, 
                  p_u, rnd1, rnd2)).block(block*dilE, block*dilE, dilE, dilE) *
               vdaggerv.return_rvdaggervr(p_d, t_sink, dirac_d, rnd2, rnd1)
                                  .block(0, block*dilE, dilE, dilE)).trace();
          }
        }} // Loops over random vectors end here! 
      }}// Loops over dirac_d and p_d end here
    }}// Loops over dirac_u and p_u end here
  }}// Loops over displacements end here

  // Using the dagger operation to get all possible random vector combinations
  // TODO: Think about imaginary correlations functions - There might be an 
  //       additional minus sign involved
  for(int displ_u = 0; displ_u < nb_dis; displ_u++){
  for(int displ_d = 0; displ_d < nb_dis; displ_d++){
    for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){
    for(int p_u = 0; p_u < nb_mom; ++p_u) {
      for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p_d = 0; p_d < nb_mom; ++p_d) {
        // TODO: A collpase of both random vectors might be better but
        //       must be done by hand because rnd2 starts from rnd1+1
        #pragma omp parallel for schedule(dynamic)
        for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
        for(int rnd2 = rnd1+1; rnd2 < nb_rnd; ++rnd2){
          Corr[p_u][p_d][dirac_d][dirac_u][displ_u]
              [displ_d][t_source][t_sink][rnd2][rnd1] = 
                   std::conj(Corr[p_u][p_d][dirac_u][dirac_d][displ_u]
                                 [displ_d][t_source][t_sink][rnd1][rnd2]); 
        }} // Loops over random vectors end here! 
      }}// Loops over dirac_d and p_d end here
    }}// Loops over dirac_u and p_u end here
  }}// Loops over displacements end here
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// build 2pt-function C2_mes for pi^+ from Corr. Equivalent to just summing
// up traces with same time difference between source and sink (all to all)
// for every dirac structure, momentum, displacement
void LapH::Correlators::build_and_write_2pt(const size_t config_i){

  std::cout << "\n\tcomputing the 2pt function:\r";
  clock_t time = clock();
  // TODO: Check if really everything is needed
  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const int max_mom_squared = global_data->get_number_of_max_mom();
  const int nb_mom = global_data->get_number_of_momenta();
  const std::vector<int> mom_squared = global_data->get_momentum_squared();
  const int p_min = 0; 
  const int p_max = nb_mom;

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // compute the norm for 4pt functions
  int norm = nb_rnd*(nb_rnd-1)*(nb_rnd-2);
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();
  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  // setting the correlation function to zero
  std::fill(C2_mes.origin(), C2_mes.origin() + C2_mes.num_elements() , 
                                                            cmplx(0.0, 0.0));
  for(int t_source = 0; t_source < Lt; ++t_source){
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    for(int p2 = 0; p2 <= max_mom_squared; p2++){
    for(int p = p_min; p < p_max; ++p){
    if(mom_squared[p] == p2){
      for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){
      for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
        for(int displ_u = 0; displ_u < nb_dis; displ_u++){
        for(int displ_d = 0; displ_d < nb_dis; displ_d++){   
          for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
          for(int rnd2 = 0; rnd2 < nb_rnd; ++rnd2){
          if(rnd1 != rnd2){
            C2_mes[p2][dirac_u][dirac_d][displ_u][displ_d]
                  [abs((t_sink - t_source - Lt) % Lt)] += 
               Corr[p][nb_mom - p - 1][dirac_u][dirac_d]
                   [displ_u][displ_d][t_source][t_sink][rnd1][rnd2];
          }}}
        }}
      }}
    }}}
  }}
  // normalization of correlation function
  double norm3 = Lt * nb_rnd * (nb_rnd - 1);
  for(auto i = C2_mes.data(); i < (C2_mes.data()+C2_mes.num_elements()); i++)
    *i /= norm3;

  // output to binary file - only diagaonal and summed momenta
  for(int dirac = 0; dirac < nb_dir; ++dirac){
    for(int p = 0; p <= max_mom_squared; p++){
      for(int displ = 0; displ < nb_dis; ++displ){
        sprintf(outfile, 
            "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
            "C2_pi+-_conf%04d.dat", 
            outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
            displ_min, displ_max, (int)config_i);
        if((fp = fopen(outfile, "wb")) == NULL)
          std::cout << "fail to open outputfile: " << outfile << std::endl;

        fwrite((double*) &(C2_mes[p][dirac][dirac][displ][displ][0]), 
                                                sizeof(double), 2 * Lt, fp);
        fclose(fp);
      }
    }
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// TODO: See if there are common parts which can be shifted to a seperate
//       function
void LapH::Correlators::build_and_write_C4_1(const size_t config_i){

  std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
  clock_t time = clock();

  // TODO: Check if really everything is needed
  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const int max_mom_squared = global_data->get_number_of_max_mom();
  const int nb_mom = global_data->get_number_of_momenta();
  const std::vector<int> mom_squared = global_data->get_momentum_squared();
  const int p_min = 0; 
  const int p_max = nb_mom;

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // compute the norm for 4pt functions
  int norm = nb_rnd*(nb_rnd-1)*(nb_rnd-2);
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();
  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  // setting the correlation function to zero
  std::fill(C4_mes.data(), C4_mes.data() + C4_mes.num_elements(), 
                                                            cmplx(0.0, 0.0));
  for(int t_source = 0; t_source < Lt; ++t_source){
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    int t_source_1 = (t_source + 1) % Lt;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){     
    for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p_u = p_min; p_u < p_max; ++p_u) {
      for(int p_d = p_min; p_d < p_max; ++p_d) {
      if(mom_squared[p_u] <= mom_squared[p_d]){
        for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
        for(int rnd2 = 0; rnd2 < nb_rnd; ++rnd2){      
        if(rnd2 != rnd1){
          for(int rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
          if((rnd3 != rnd2) && (rnd3 != rnd1)){
            for(int rnd4 = 0; rnd4 < nb_rnd; ++rnd4){
            if((rnd4 != rnd1) && (rnd4 != rnd2) && (rnd4 != rnd3)){
              C4_mes[p_u][p_d][dirac_u][dirac_d]
                    [abs((t_sink - t_source - Lt) % Lt)] +=
                (Corr[p_u][nb_mom - p_d - 1][dirac_u][dirac_d]
                     [0][0][t_source_1][t_sink_1][rnd1][rnd3]) *
                (Corr[nb_mom - p_u - 1][p_d][dirac_u][dirac_d]
                     [0][0][t_source][t_sink][rnd2][rnd4]);
            }}// loop rnd4 and if
          }}// loop rnd3 and if
        }}}// loops rnd1 and rnd 2 and if
      }}}// loops momenta and if
    }}// loops dirac
  }}// loops t_sink and t_source
  // Normalization of 4pt-function
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to binary file
  // see output to binary file for C2. 
  // write into folders with suffix "_unsuppressed". These only include
  // correlators of the diagonal matrix elements of the GEVP for which
  // the three-momentum remains unchanged for both quarks. Because the
  // quarks have to be back-to-back, for the offdiagonal elements this
  // cannot occur. The suppression can be interpreted as Zweig-suppressed
  // gluon exchange
  for(int dirac = 0; dirac < nb_dir; ++dirac){
    for(int p = 0; p <= max_mom_squared; p++){
      sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d_unsuppressed/"
          "C4_1_conf%04d.dat", 
          outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
          displ_min, displ_max, (int)config_i);
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
  for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){
    for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p1 = 0; p1 <= max_mom_squared; p1++){
        for(int p2 = p1; p2 <= max_mom_squared; p2++){

          sprintf(outfile, 
             "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
             "C4_1_conf%04d.dat", 
             outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
             p1, p2, displ_min, displ_max,(int) config_i);
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
  time = clock() - time;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::build_and_write_C4_2(const size_t config_i){

  std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
  clock_t time = clock();

  // TODO: Check if really everything is needed
  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const int max_mom_squared = global_data->get_number_of_max_mom();
  const int nb_mom = global_data->get_number_of_momenta();
  const std::vector<int> mom_squared = global_data->get_momentum_squared();
  const int p_min = 0; 
  const int p_max = nb_mom;

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // compute the norm for 4pt functions
  int norm = nb_rnd*(nb_rnd-1)*(nb_rnd-2);
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  const int dirac_min = global_data->get_dirac_min();
  const int dirac_max = global_data->get_dirac_max();
  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();
  // setting the correlation function to zero
  std::fill(C4_mes.data(), C4_mes.data() + C4_mes.num_elements(), 
                                                             cmplx(0.0, 0.0));
  for(int t_source = 0; t_source < Lt; ++t_source){
  for(int t_sink = 0; t_sink < Lt - 1; ++t_sink){
    int t_source_1 = (t_source + 1) % Lt;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){     
    for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p_u = p_min; p_u < p_max; ++p_u) {
      for(int p_d = p_min; p_d < p_max; ++p_d) {
      if(mom_squared[p_u] <= mom_squared[p_d]){
        for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
        for(int rnd2 = 0; rnd2 < nb_rnd; ++rnd2){      
        if(rnd2 != rnd1){
          for(int rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
          if((rnd3 != rnd2) && (rnd3 != rnd1)){
            for(int rnd4 = 0; rnd4 < nb_rnd; ++rnd4){
            if((rnd4 != rnd1) && (rnd4 != rnd2) && (rnd4 != rnd3)){
              C4_mes[p_u][p_d][dirac_u][dirac_d]
                    [abs((t_sink - t_source - Lt) % Lt)] +=
                (Corr[p_u][nb_mom - p_d - 1][dirac_u][dirac_d]
                     [0][0][t_source_1][t_sink][rnd1][rnd3]) *
                (Corr[nb_mom - p_u - 1][p_d][dirac_u][dirac_d]
                     [0][0][t_source][t_sink_1][rnd2][rnd4]);
            }}// loop rnd4
          }}// loop rnd3
        }}}// loops rnd2 and rnd1
      }}}// loops momenta
    }}// loops dirac
  }}// loops t_source and t_sink
  // Normalization of 4pt-function. 
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to binary file
  // see output to binary file for C2. 
  // write into folders with suffix "_unsuppressed". These only include
  // correlators of the diagonal matrix elements of the GEVP for which
  // the three-momentum remains unchanged for both quarks. Because the
  // quarks have to be back-to-back, for the offdiagonal elements this
  // cannot occur. The suppression can be interpreted as Zweig-suppressed
  // gluon exchange
  for(int dirac = 0; dirac < nb_dir; ++dirac){
    for(int p = 0; p <= max_mom_squared; p++){
      sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d_unsuppressed/"
          "C4_2_conf%04d.dat", 
          outpath.c_str(), dirac_ind.at(dirac), dirac_ind.at(dirac), p, p, 
          displ_min, displ_max, (int)config_i);
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

  for(int dirac_u = 0; dirac_u < nb_dir; ++dirac_u){
    for(int dirac_d = 0; dirac_d < nb_dir; ++dirac_d){
      for(int p1 = 0; p1 <= max_mom_squared; p1++){
        for(int p2 = p1; p2 <= max_mom_squared; p2++){

          sprintf(outfile, 
             "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
             "C4_2_conf%04d.dat", 
             outpath.c_str(), dirac_ind.at(dirac_u), dirac_ind.at(dirac_d), 
             p1, p2, displ_min, displ_max, (int)config_i);
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
  time = clock() - time;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
}
