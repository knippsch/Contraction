#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::Correlators::build_Corr(){

  const size_t Lt = global_data->get_Lt();
  const vec_pdg_C2 op_C2 = global_data->get_op_C2();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int dilT = quarks[0].number_of_dilution_T;
  const indexlist_2 rnd_vec_index = global_data->get_rnd_vec_C2();
  // TODO: must be changed in GlobalData {
  // TODO: }

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
    //#pragma omp parallel
    //#pragma omp single
    //{
      for(const auto& op : op_C2){
      for(const auto& i : op.index){
    
        // TODO: A collpase of both random vectors might be better but
        //       must be done by hand because rnd2 starts from rnd1+1
        for(const auto& rnd_it : rnd_vec_index) {
          // build all 2pt traces leading to C2_mes
          // Corr = tr(D_d^-1(t_sink) Gamma D_u^-1(t_source) Gamma)
          // TODO: Just a workaround
    
          //#pragma omp task shared(rnd_it, i)
          compute_meson_small_traces(i.second, basic.get_operator
              (t_source, t_sink/dilT, 1, i.first, rnd_it.first, rnd_it.second),
            //TODO: shouldn't that be op_Corr[i.second].id_rVdaggerVr?
            vdaggerv.return_rvdaggervr(i.second, t_sink, rnd_it.second, rnd_it.first),
            Corr[i.first][i.second][t_source][t_sink][rnd_it.first][rnd_it.second]);

        } // Loop over random vectors ends here! 
      }}//Loops over all Quantum numbers 
    
      // Using the dagger operation to get all possible random vector combinations
      // TODO: Think about imaginary correlations functions - There might be an 
      //       additional minus sign involved
    //} // end parallel region
    
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
  // TODO: Check if really everything is needed
  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const vec_pdg_C2 op_C2 = global_data->get_op_C2();
  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  const size_t nb_op = global_data->get_number_of_operators();
  const indexlist_2 rnd_vec_index = global_data->get_rnd_vec_C2();

  // compute the norm for 2pt functions
  int norm = rnd_vec_index.size();
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  // setting the correlation function to zero
  std::fill(C2_mes.origin(), C2_mes.origin() + C2_mes.num_elements() , 
                                                            cmplx(0.0, 0.0));
  for(int t_source = 0; t_source < Lt; ++t_source){
  for(int t_sink = 0; t_sink < Lt; ++t_sink){

    for(const auto& op : op_C2) {
    for(const auto& i : op.index) {
      for(const auto& rnd : rnd_vec_index) {
//        C2_mes[op.p_sq][op.dg_so][op.dg_si]
        C2_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] += 
           Corr[i.first][i.second][t_source][t_sink][rnd.first][rnd.second];
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

  sprintf(outfile, 
      "%s/C2_pi+-_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_2pt(outfile, C2_mes);

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

  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();
  const int nb_mom = global_data->get_number_of_momenta();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();

  // compute the norm for 4pt functions
  int norm = rnd_vec_index.size();
  std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
  const double norm1 = Lt * norm;

  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

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

    for(const auto& op : op_C4){
    for(const auto& i : op.index){
      for(const auto& rnd : rnd_vec_index) {
        C4_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] +=
          (Corr[i[0]][i[2]][t_source_1][t_sink_1][rnd[0]][rnd[2]]) *
          (Corr[i[1]][i[3]][t_source][t_sink][rnd[1]][rnd[3]]);
      } // loop over random vectors
    }}//loops operators
  }}// loops t_sink and t_source
  // Normalization of 4pt-function
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // to build a GEVP, the correlators are written into a seperate folder
  // for every dirac structure, momentum, (entry of the GEVP matrix).
  // displacement is not supported at the moment
  for(const auto& op : op_C4){
    size_t gam_so = op_Corr[op.index.front()[0]].gamma[0];
    size_t gam_si = op_Corr[op.index.front()[2]].gamma[0];
    size_t mom_so = op_Corr[op.index.front()[0]].p3[0] * op_Corr[op.index.front()[0]].p3[0] + 
                  op_Corr[op.index.front()[0]].p3[1] * op_Corr[op.index.front()[0]].p3[1] + 
                  op_Corr[op.index.front()[0]].p3[2] * op_Corr[op.index.front()[0]].p3[2];
    size_t mom_si = op_Corr[op.index.front()[2]].p3[0] * op_Corr[op.index.front()[2]].p3[0] + 
                  op_Corr[op.index.front()[2]].p3[1] * op_Corr[op.index.front()[2]].p3[1] + 
                  op_Corr[op.index.front()[2]].p3[2] * op_Corr[op.index.front()[2]].p3[2];
    size_t dis_so = op_Corr[op.index.front()[0]].dis3[0] * op_Corr[op.index.front()[0]].dis3[0] + 
                  op_Corr[op.index.front()[0]].dis3[1] * op_Corr[op.index.front()[0]].dis3[1] + 
                  op_Corr[op.index.front()[0]].dis3[2] * op_Corr[op.index.front()[0]].dis3[2];
    size_t dis_si = op_Corr[op.index.front()[2]].dis3[0] * op_Corr[op.index.front()[2]].dis3[0] + 
                  op_Corr[op.index.front()[2]].dis3[1] * op_Corr[op.index.front()[2]].dis3[1] + 
                  op_Corr[op.index.front()[2]].dis3[2] * op_Corr[op.index.front()[2]].dis3[2];

    sprintf(outfile, 
        "%s/dirac_%02ld_%02ld_p_%01ld_%01ld_displ_%01ld_%01ld/"
        "C4_1_conf%04d.dat", 
        outpath.c_str(), gam_so, gam_si, mom_so, mom_si, 
        dis_so, dis_si, (int)config_i);
    std::cout << outfile << std::endl;
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile: " << outfile << std::endl;

    fwrite((double*) &(C4_mes[op.id][0]), 
                                            sizeof(double), 2 * Lt, fp);
    fclose(fp);


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

  // global variables from input file needed here
  const int Lt = global_data->get_Lt();

  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();
  const int nb_mom = global_data->get_number_of_momenta();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();

 // compute the norm for 4pt functions
  int norm = rnd_vec_index.size();
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

    for(const auto& op : op_C4){
    for(const auto& i : op.index){
      for(const auto& rnd : rnd_vec_index) {
        C4_mes[op.id][abs((t_sink - t_source - Lt) % Lt)] +=
          (Corr[i[0]][i[2]][t_source_1][t_sink][rnd[0]][rnd[2]]) *
          (Corr[i[1]][i[3]][t_source][t_sink_1][rnd[1]][rnd[3]]);
      } // loop over random vectors
    }}//loops operators
  }}// loops t_source and t_sink

  // Normalization of 4pt-function. 
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // to build a GEVP, the correlators are written into a seperate folder
  // for every dirac structure, momentum, (entry of the GEVP matrix).
  // displacement is not supported at the moment
  for(const auto& op : op_C4){
    size_t gam_so = op_Corr[op.index.front()[0]].gamma[0];
    size_t gam_si = op_Corr[op.index.front()[2]].gamma[0];
    size_t mom_so = op_Corr[op.index.front()[0]].p3[0] * op_Corr[op.index.front()[0]].p3[0] + 
                  op_Corr[op.index.front()[0]].p3[1] * op_Corr[op.index.front()[0]].p3[1] + 
                  op_Corr[op.index.front()[0]].p3[2] * op_Corr[op.index.front()[0]].p3[2];
    size_t mom_si = op_Corr[op.index.front()[2]].p3[0] * op_Corr[op.index.front()[2]].p3[0] + 
                  op_Corr[op.index.front()[2]].p3[1] * op_Corr[op.index.front()[2]].p3[1] + 
                  op_Corr[op.index.front()[2]].p3[2] * op_Corr[op.index.front()[2]].p3[2];
    size_t dis_so = op_Corr[op.index.front()[0]].dis3[0] * op_Corr[op.index.front()[0]].dis3[0] + 
                  op_Corr[op.index.front()[0]].dis3[1] * op_Corr[op.index.front()[0]].dis3[1] + 
                  op_Corr[op.index.front()[0]].dis3[2] * op_Corr[op.index.front()[0]].dis3[2];
    size_t dis_si = op_Corr[op.index.front()[2]].dis3[0] * op_Corr[op.index.front()[2]].dis3[0] + 
                  op_Corr[op.index.front()[2]].dis3[1] * op_Corr[op.index.front()[2]].dis3[1] + 
                  op_Corr[op.index.front()[2]].dis3[2] * op_Corr[op.index.front()[2]].dis3[2];

    sprintf(outfile, 
        "%s/dirac_%02ld_%02ld_p_%01ld_%01ld_displ_%01ld_%01ld/"
        "C4_2_conf%04d.dat", 
        outpath.c_str(), gam_so, gam_si, mom_so, mom_si, 
        dis_so, dis_si, (int)config_i);
    std::cout << outfile << std::endl;
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile: " << outfile << std::endl;

    fwrite((double*) &(C4_mes[op.id][0]), 
                                            sizeof(double), 2 * Lt, fp);
    fclose(fp);


  }

  time = clock() - time;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
}

