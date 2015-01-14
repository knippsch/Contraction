#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C4_3() ?
void LapH::Correlators::compute_meson_4pt_cross_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();
  // TODO: must be changed in GlobalData {
  // TODO: }

  std::cout << "\n\tcomputing the traces of 2 pi_+/-:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of 2 pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){
      const int t_source_1 = (t_source + 1) % Lt;

      if(t_source != 0){
        if(t_source%2 == 0){
          X.swap(1, 0);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
        }
        else{
          X.swap(0, 1);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 0);
        }
      }
      else{
        X.construct(basic, vdaggerv, 0, t_source,   t_sink, 0);
        X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
      }
      if(t_source == t_sink)
        continue;
    
      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
        #pragma omp task shared(op)
        for(const auto& i : op.index){
        // complete diagramm. combine X and Y to four-trace
        // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
        //             D_d^-1(t_sink       | t_source + 1) Gamma 
        //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
        //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            if(t_source%2 == 0)

              priv_C4 += (X(0, i[2], i[1], rnd_it[2], rnd_it[1], rnd_it[3]) *
                          X(1, i[3], i[0], rnd_it[3], rnd_it[0], rnd_it[2])).trace();
            else
              priv_C4 += std::conj(
                         (X(0, i[2], i[1], rnd_it[2], rnd_it[1], rnd_it[3]) *
                          X(1, i[3], i[0], rnd_it[3], rnd_it[0], rnd_it[2])).trace());
          }
          #pragma omp critical
          {
            C4_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the traces of 2 pi_+/-: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void LapH::Correlators::write_C4_3(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();
  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();

  // normalisation
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to binary file
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
        "C4_3_conf%04d.dat", 
        outpath.c_str(), gam_so, gam_si, mom_so, mom_si, 
        dis_so, dis_si, (int)config_i);
    std::cout << outfile << std::endl;
    if((fp = fopen(outfile, "wb")) == NULL)
      std::cout << "fail to open outputfile: " << outfile << std::endl;

    fwrite((double*) &(C4_mes[op.id][0]), sizeof(double), 2 * Lt, fp);
    fclose(fp);
  }
}

