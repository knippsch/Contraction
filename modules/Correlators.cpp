#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
LapH::Correlators::Correlators() : basic(), peram(), rnd_vec(), vdaggerv(),
                                   C4_mes(), C2_mes(), Corr()  {

  const vec_index_2pt op_C2 = global_data->get_lookup_2pt_trace();
  const size_t nb_op_2pt = op_C2.size();
  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const size_t nb_op_4pt = op_C4.size();
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  const size_t nb_op = op_Corr.size();
//  const size_t nb_op = global_data->get_number_of_operators();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  // TODO: }

  rnd_vec.resize(nb_rnd, LapH::RandomVector(Lt*nb_ev*4));

  //TODO: size of C4_mes and C2_mes must be replaced by size of corresponding
  //operator lists. Momentary values are upper limit
  C4_mes.resize(boost::extents[nb_op_4pt][Lt]);
  C2_mes.resize(boost::extents[nb_op_2pt][Lt]);
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

  basic.init_operator('b', vdaggerv, peram);

  // computing the meson correlator which can be used to compute all small
  // trace combinations for 2pt and 4pt functions
  build_Corr();

  // computing the meson 4pt big cross trace
  // TODO: if condition that at least four random vectos are needed
//  compute_meson_4pt_cross_trace(X);
//
//  write_C4_3(config_i);
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
