#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();


LapH::Correlators::Correlators() : basic(), peram(), rnd_vec(), vdaggerv(),
                                   C4_mes(), C2_mes(), Corr()  {

  const size_t Lt = global_data->get_Lt();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  // TODO: must be changed in GlobalData {
  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

  rnd_vec.resize(nb_rnd, LapH::RandomVector(Lt*nb_ev*4));

  C4_mes.resize(boost::extents[nb_mom][nb_mom][nb_dir][nb_dir][Lt]);
  C2_mes.resize(boost::extents[nb_mom][nb_dir][nb_dir][nb_dis][nb_dis][Lt]);
  Corr.resize(boost::extents[nb_mom][nb_mom][nb_dir][nb_dir][nb_dis]
                                    [nb_dis][Lt][Lt][nb_rnd][nb_rnd]);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::compute_meson_corr(const int t_source, 
                                           const int t_sink){

  const size_t Lt = global_data->get_Lt();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t dilT = quarks[0].number_of_dilution_T;
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
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
        for(int rnd2 = rnd1+1; rnd2 < nb_rnd; ++rnd2){
          // build all 2pt traces leading to C2_mes
          // Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)
          // TODO: Just a workaround
          std::array<double, 4> bla = {{1., 1., -1., -1.}};
          for(size_t block = 0; block < 4; block++){
            // TODO: dilution scheme in time should be choosable
            const size_t so = (t_source/dilT)*4*dilE + block*dilE;
            Corr[p_u][p_d][dirac_u][dirac_d][displ_u][displ_d]
              [t_source][t_sink][rnd1][rnd2] += bla[block] *
              ((basic.get_operator(0, dirac_u, p_u, rnd1, rnd2))
                                   .block(so, so, dilE, dilE) *
               vdaggerv.return_rvdaggervr(p_d, t_source, dirac_d, rnd2, rnd1)
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
        #pragma omp parallel for collapse(1) schedule(dynamic)
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
void LapH::Correlators::compute_meson_4pt_cross(LapH::CrossOperator& X,
                                                const int t_source, 
                                                const int t_sink){
  const int Lt = global_data->get_Lt();
  const int t_source_1 = (t_source + 1) % Lt;
  const int t_sink_1 = (t_sink + 1) % Lt;
  const size_t nb_mom = global_data->get_number_of_momenta();
  const int max_mom_squared = global_data->get_number_of_max_mom();
  const std::vector<int> mom_squared = global_data->get_momentum_squared();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  // TODO: must be changed in GlobalData {
  int displ_min = global_data->get_displ_min();
  int displ_max = global_data->get_displ_max();
  const size_t nb_dis = displ_max - displ_min + 1;
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

//  if(t_source != 0){
//    X.swap(0, 1);
//    X.construct(basic, vdaggerv, 1, t_source, 1);
//  }
//  else{
//    X.construct(basic, vdaggerv, 0, t_source, 0);
//    X.construct(basic, vdaggerv, 1, t_source, 1);
//  }
//  if(t_source == t_sink)
//    return;

  if(t_source == t_sink)
    return;
  X.construct(basic, vdaggerv, 0, t_source, 0);
  X.construct(basic, vdaggerv, 1, t_source, 1);

  for(size_t dirac_1 = 0; dirac_1 < nb_dir; ++dirac_1){     
    for(size_t p = 0; p <= max_mom_squared; p++){
    for(size_t p_u = 0; p_u < nb_mom; ++p_u) {
    if(mom_squared[p_u] == p){
        for(size_t dirac_2 = 0; dirac_2 < nb_dir; ++dirac_2){
        for(size_t p_d = 0; p_d < nb_mom; ++p_d) {
        if(mom_squared[p_d] == p){
          // complete diagramm. combine X and Y to four-trace
          // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
          //             D_d^-1(t_sink       | t_source + 1) Gamma 
          //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
          //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          #pragma omp parallel
          {
            cmplx priv_C4(0.0,0.0);
            #pragma omp for collapse(2) schedule(dynamic)
            for(size_t rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
            for(size_t rnd2 = 0; rnd2 < nb_rnd; ++rnd2){      
            if(rnd2 != rnd1){
            for(size_t rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
            if((rnd3 != rnd2) && (rnd3 != rnd1)){
            for(size_t rnd4 = 0; rnd4 < nb_rnd; ++rnd4){
            if((rnd4 != rnd1) && (rnd4 != rnd2) && (rnd4 != rnd3)){
// TODO: think about dirac structure for off-diagonal elements
//              if(t_source%2 == 0)
                priv_C4 += (X(0, p_d, p_u, dirac_1, dirac_2, rnd3, rnd2, rnd4) *
                            X(1, nb_mom - p_d - 1, nb_mom - p_u - 1,
                              dirac_1, dirac_2, rnd4, rnd1, rnd3)).trace();
//              else
//                priv_C4 += std::conj(
//                           (X(0, p_d, p_u, dirac_1, dirac_2, rnd3, rnd2, rnd4) * 
//                            X(1, nb_mom - p_d - 1, nb_mom - p_u - 1,
//                              dirac_1, dirac_2, rnd4, rnd1, rnd3)).trace());
            }}}}}}}
            #pragma omp critical
            {
              C4_mes[p][p][dirac_1][dirac_2]
                  [abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
            }
          }
        }}// loop and if condition p_d
      }// loop dirac_2
    }}}// loop and if conditions p_u
  }// loop dirac_1
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::build_everything(const size_t config_i){

  // initialising the big arrays
  set_corr(config_i);
  // setting the correlation functions to zero
  std::fill(Corr.data(), Corr.data()+Corr.num_elements(), cmplx(.0,.0));
  std::fill(C4_mes.data(), C4_mes.data()+C4_mes.num_elements(), cmplx(.0,.0));

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

  // compute the norm for 4pt functoins
  int norm = 0;
  for(int rnd1 = 0; rnd1 < nb_rnd; ++rnd1)
    for(int rnd2 = 0; rnd2 < nb_rnd; ++rnd2)      
      if(rnd2 != rnd1)
        for(int rnd3 = 0; rnd3 < nb_rnd; ++rnd3)
          if((rnd3 != rnd2) && (rnd3 != rnd1))
            for(int rnd4 = 0; rnd4 < nb_rnd; ++rnd4)
              if((rnd4 != rnd1) && (rnd4 != rnd2) && (rnd4 != rnd3))
                norm++;
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;


  // memory for intermediate matrices when building C4_3 (save multiplications)
  LapH::CrossOperator X(2);

  std::cout << "\n\tcomputing the traces of pi_+/-:\r";
  clock_t time = clock();
  // calculate all two-operator traces of the form tr(u \Gamma \bar{d}) and 
  // build all combinations of momenta, dirac_structures and displacements as
  // specified in infile






  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    // this is a way to reuse the already computed operators from t_source_1
    // initialize contraction[rnd_i] = perambulator * basicoperator = D_u^-1
    // choose 'i' for interlace or 'b' for block time dilution scheme
    // TODO: get that from input file
    if(t_sink != 0){
      basic.swap_operators();
      basic.init_operator(1, t_sink_1, 'b', 0, vdaggerv, peram);
    }
    else {
      basic.init_operator(0, t_sink,   'b', 0, vdaggerv, peram);
      basic.init_operator(1, t_sink_1, 'b', 0, vdaggerv, peram);
    }
    for(int t_source = 0; t_source < Lt; ++t_source){
      // computing the meson correlator which can be used to compute all small
      // trace combinations for 2pt and 4pt functions
      compute_meson_corr(t_source, t_sink);
      // computing the meson 4pt big cross trace
      // TODO: if condition that at least four random vectos are needed
      compute_meson_4pt_cross(X, t_source, t_sink);
    }
  }// Loops over time end here










  // *************************************************************************
  // FOUR PT CONTRACTION 3 ***************************************************
  // *************************************************************************
  // Normalization of 4pt-function
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;
  // output to binary file
  for(size_t dirac_1 = 0; dirac_1 < nb_dir; ++dirac_1){     
  for(size_t dirac_2 = 0; dirac_2 < nb_dir; ++dirac_2){
    for(size_t p = 0; p <= max_mom_squared; p++){
      sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
          "C4_3_conf%04d.dat", 
          outpath.c_str(), dirac_ind.at(dirac_1), dirac_ind.at(dirac_2), 
          (int)p, (int)p, 0, 0, (int)config_i);
      if((fp = fopen(outfile, "wb")) == NULL)
        std::cout << "fail to open outputfile" << std::endl;
      fwrite((double*) &(C4_mes[p][p][dirac_1][dirac_2][0]), 
             sizeof(double), 2 * Lt, fp);
      fclose(fp);
    }// loop p
  }}// loops dirac_1 dirac_2

  ////////////////////////////////////////////////////////////////////////////
  //                          TWO POINT FUNCTION                            //
  ////////////////////////////////////////////////////////////////////////////
  // build 2pt-function C2_mes for pi^+ from Corr. Equivalent to just summing
  // up traces with same time difference between source and sink (all to all)
  // for every dirac structure, momentum, displacement
  std::cout << "\tcomputing the pi_+/-\n" << std::endl;
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


  ////////////////////////////////////////////////////////////////////////////
  //                         FOUR POINT FUNCTION                            //
  ////////////////////////////////////////////////////////////////////////////
  // *************************************************************************
  // FOUR PT CONTRACTION 1 ***************************************************
  // *************************************************************************
  std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
  time = clock();
  // displacement not supported for 4pt functions atm
  displ_min = 0;
  displ_max = 0;
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


  // *************************************************************************
  // FOUR PT CONTRACTION 2 ***************************************************
  // *************************************************************************
  std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
  time = clock();
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
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::read_rnd_vectors_from_file (const int config_i) {

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
