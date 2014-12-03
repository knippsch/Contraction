#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
LapH::Correlators::Correlators() : basic(), peram(), rnd_vec(), vdaggerv(),
                                   C4_mes(), C2_mes(), Corr()  {

  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const size_t nb_op = global_data->get_number_of_operators();
  const size_t nb_dg = global_data->get_number_of_displ_gamma();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  // TODO: }

  rnd_vec.resize(nb_rnd, LapH::RandomVector(Lt*nb_ev*4));

  C4_mes.resize(boost::extents[nb_mom_sq][nb_mom_sq][nb_dg][nb_dg][Lt]);
  C2_mes.resize(boost::extents[nb_mom_sq][nb_dg][nb_dg][Lt]);
  Corr.resize(boost::extents[nb_op][nb_op][Lt][Lt][nb_rnd][nb_rnd]);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::compute_correlators(const size_t config_i){

  // initialising the big arrays
  set_corr(config_i);
  // setting the correlation functions to zero
  std::fill(Corr.data(), Corr.data()+Corr.num_elements(), cmplx(.0,.0));
  std::fill(C4_mes.data(), C4_mes.data()+C4_mes.num_elements(), cmplx(.0,.0));

  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  // memory for intermediate matrices when building C4_3 (save multiplications)
  LapH::CrossOperator X(2);

  basic.init_operator('b', 0, vdaggerv, peram);

  std::cout << "\n\tcomputing the traces of pi_+/-:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){
      // computing the meson correlator which can be used to compute all small
      // trace combinations for 2pt and 4pt functions
      compute_meson_small_traces(t_source, t_sink);
      // computing the meson 4pt big cross trace
      // TODO: if condition that at least four random vectos are needed
      compute_meson_4pt_cross_trace(X, t_source, t_sink);
    }
  }// Loops over time end here
  time = clock() - time;
  std::cout << "\n\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

  write_C4_3(config_i);
  build_and_write_2pt(config_i);
  build_and_write_C4_1(config_i);
  build_and_write_C4_2(config_i);

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
