#include "Correlators.h"

//TODO outfile can be hidden in the export_meson_2pt() (and rename that thing)

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::build_Corr(){

  const size_t Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int dilT = quarks[0].number_of_dilution_T;

  const vec_index_2pt op_C2 = global_data->get_lookup_2pt_trace();
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  const indexlist_2 rnd_vec_index = global_data->get_rnd_vec_2pt();
  // TODO: must be changed in GlobalData {
  // TODO: }     Was ist das? - Markus

  std::cout << "\n\tcomputing the traces of pi_+/-:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){

    // CJ: OMP cannot do the parallelization over an auto loop,
    //     implementing a workaround by hand
    //     not used at the moment, increases runtime by a factor of 2
    #pragma omp parallel
    #pragma omp single
    {
      for(const auto& op : op_C2){

        size_t id_Q2 = op.index_Q2;
        size_t id_Corr = op.index_Corr;
      #pragma omp task shared(op)
    
        // TODO: A collpase of both random vectors might be better but
        //       must be done by hand because rnd2 starts from rnd1+1
        for(const auto& rnd_it : rnd_vec_index) {
          // build all 2pt traces leading to C2_mes
          // Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)
          // TODO: Just a workaround
    
          compute_meson_small_traces(id_Corr, basic.get_operator
            (t_source, t_sink/dilT, 1, id_Q2, rnd_it.first, rnd_it.second),
            vdaggerv.return_rvdaggervr(op_Corr[id_Corr].id_rvdvr, t_sink, 
                rnd_it.second, rnd_it.first),
            Corr[id_Q2][id_Corr][t_source][t_sink]
                [rnd_it.first][rnd_it.second]);

        } // Loop over random vectors ends here! 
//      #pragma omp taskwait
      }//Loops over all Quantum numbers 
    
      // Using the dagger operation to get all possible random vector combinations
      // TODO: Think about imaginary correlations functions - There might be an 
      //       additional minus sign involved
    } // end parallel region
    
    }// Loops over t_source
  }// Loops over t_sink

  std::cout << "\tcomputing the traces of pi_+/-: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//calculate tr(2pt function). Only diagonal blocks are needed.
// Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)

void LapH::Correlators::compute_meson_small_traces(const size_t id_si, 
                                             const Eigen::MatrixXcd& Q2,
                                             const Eigen::MatrixXcd& rVdaggerVr,
                                             cmplx& Corr) {

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;

  for(size_t block = 0; block < 4; block++){

    cmplx value = 1;
    basic.value_dirac(id_si, block, value);

    Corr += value * (Q2.block(block*dilE, block*dilE, dilE, dilE) *
            rVdaggerVr.block(0, (basic.order_dirac(id_si, block)*dilE), 
            dilE, dilE)).trace();
    }
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

  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const vec_index_2pt op_C2 = global_data->get_lookup_2pt_trace();
  const vec_index_IO_1 op_C2_IO = global_data->get_lookup_2pt_IO();
  const indexlist_2 rnd_vec_index = global_data->get_rnd_vec_2pt();

  // compute the norm for 2pt functions
  int norm = rnd_vec_index.size();
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  // setting the correlation function to zero
  std::fill(C2_mes.origin(), C2_mes.origin() + C2_mes.num_elements() , 
                                                            cmplx(0.0, 0.0));
  for(int t_source = 0; t_source < Lt; ++t_source){
  for(int t_sink = 0; t_sink < Lt; ++t_sink){

    for(const auto& op : op_C2_IO) {
    for(const auto& i : op.index_pt){
      size_t id_Q2 = op_C2[i].index_Q2;
      size_t id_Corr = op_C2[i].index_Corr;

      for(const auto& rnd : rnd_vec_index) {
        C2_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] = 
           Corr[id_Q2][id_Corr][t_source][t_sink][rnd.first][rnd.second];
      } //Loop over random vectors
    }}//Loops over all Quantum numbers
  }}//Loops over time

  // normalization of correlation function
  for(auto i = C2_mes.data(); i < (C2_mes.data()+C2_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile     - filename
  // run_id      - first message with runinfo specific for each run of the 
  //               program
  // attributes  - vector of tags containing quantum numbers for each correlator
  // correlators - vector of correlators

  sprintf(outfile, "%s/C2_pi+-_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C2_IO, "C2+", C2_mes);

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

  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

//  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const vec_index_2pt op_C2 = global_data->get_lookup_2pt_trace();
  const vec_index_IO_2 op_C4_IO = global_data->get_lookup_4pt_1_IO();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();

  // compute the norm for 4pt functions
  int norm = rnd_vec_index.size();
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

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

    for(const auto& op : op_C4_IO){
    for(const auto& i : op.index_pt){
      size_t id_Q2_0   = op_C2[i.first].index_Q2;
      size_t id_Corr_0 = op_C2[i.first].index_Corr;
      size_t id_Q2_1   = op_C2[i.second].index_Q2;
      size_t id_Corr_1 = op_C2[i.second].index_Corr;

      for(const auto& rnd : rnd_vec_index) {
        C4_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] +=
          (Corr[id_Q2_0][id_Corr_0][t_source_1][t_sink_1][rnd[0]][rnd[2]]) *
          (Corr[id_Q2_1][id_Corr_1][t_source][t_sink][rnd[1]][rnd[3]]);
      } // loop over random vectors
    }}//loops operators
  }}// loops t_sink and t_source
  // Normalization of 4pt-function
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

// output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4_1_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I2+_1", C4_mes);

  time = clock() - time;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void LapH::Correlators::build_and_write_C4_2(const size_t config_i){

  std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
  clock_t time = clock();

  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

//  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const vec_index_2pt op_C2 = global_data->get_lookup_2pt_trace();
  const vec_index_IO_2 op_C4_IO = global_data->get_lookup_4pt_2_IO();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();

 // compute the norm for 4pt functions
  int norm = rnd_vec_index.size();
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

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

    for(const auto& op : op_C4_IO){
    for(const auto& i : op.index_pt){
      size_t id_Q2_0   = op_C2[i.first].index_Q2;
      size_t id_Corr_0 = op_C2[i.first].index_Corr;
      size_t id_Q2_1   = op_C2[i.second].index_Q2;
      size_t id_Corr_1 = op_C2[i.second].index_Corr;

      for(const auto& rnd : rnd_vec_index) {
        C4_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] +=
          (Corr[id_Q2_0][id_Corr_0][t_source_1][t_sink][rnd[0]][rnd[2]]) *
          (Corr[id_Q2_1][id_Corr_1][t_source][t_sink_1][rnd[1]][rnd[3]]);
      } // loop over random vectors
    }}//loops operators
  }}// loops t_source and t_sink

  // Normalization of 4pt-function. 
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4_2_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I2+_2", C4_mes);

  time = clock() - time;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
}

