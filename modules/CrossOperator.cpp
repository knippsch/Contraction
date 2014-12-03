#include "CrossOperator.h"  

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

LapH::CrossOperator::CrossOperator(const size_t number) : X(number) {

  const size_t nb_mom = global_data->get_number_of_momenta();
  const size_t nb_op = nb_mom;                                                //!!!!
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // TODO: must be changed in GlobalData {
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

  for(auto& xx : X){
    xx.resize(boost::extents[nb_op][nb_op][nb_rnd][nb_rnd][nb_rnd]);

    std::fill(xx.data(), xx.data() + xx.num_elements(), 
                Eigen::MatrixXcd(4 * dilE, 4 * dilE));
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::construct(const BasicOperator& basic, 
                                    const VdaggerV& vdaggerv, const size_t nb,
                                    const int t_source, const int t_sink,
                                    const size_t type){

  const int Lt = global_data->get_Lt();
  const size_t nb_dg = global_data->get_number_of_displ_gamma();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t dilT = quarks[0].number_of_dilution_T;
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();
  // TODO: must be changed in GlobalData {
  // TODO: }

  size_t tu, td, t2;
  if(type == 0){
    tu = (t_sink/dilT);
    if (tu == (((t_sink+1)%Lt)/dilT))
      td = 1;
    else
      td = 2;
    t2 = (t_sink + 1)%Lt;
  }
  else{
    t2 = (t_sink + 1)%Lt;
    tu = (t2/dilT);
    if (tu == t_sink/dilT)
      td = 1;
    else
      td = 0;
    t2 = t_sink;
  }

#pragma omp parallel
{     
  Eigen::MatrixXcd rvdaggervr = Eigen::MatrixXcd::Zero(dilE, 4*dilE);

  for(const auto& op : op_C4){
  for(const auto& i : op.index){
    size_t id_so = i[2+nb];
    size_t id_si = i[1-nb];

    #pragma omp for collapse(2) schedule(dynamic)
    for(size_t rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
    for(size_t rnd2 = 0; rnd2 < nb_rnd; ++rnd2){
    if(rnd2 != rnd1){
    for(size_t rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
    if((rnd3 != rnd1) && (rnd3 != rnd2)){

      basic.mult_dirac(vdaggerv.return_rvdaggervr(id_si/nb_dg, t2, op.dg_si, 
          rnd2, rnd3), rvdaggervr, id_si);
      for(size_t block = 0; block < 4; block++){

        X[nb][id_so][id_si][rnd1][rnd2][rnd3]
                                .block(0, block*dilE, 4*dilE, dilE) = 
          basic.get_operator(t_source, tu, td, id_so, rnd1, rnd2)
                                .block(0, block*dilE, 4*dilE, dilE) * 
          rvdaggervr.block(0, block*dilE, dilE, dilE);

      }// loop block ends here
    }}}}}// loops random vectors
  }}// loops operators
}

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::swap(const size_t nb1, const size_t nb2){
  
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  // TODO: }
  // TODO: Think about for each loop

  // can one collapse autoloops?
  for(auto& op_so : op_Corr){
  for(auto& op_si : op_Corr){
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
    for(size_t rnd2 = 0; rnd2 < nb_rnd; ++rnd2){
    if(rnd2 != rnd1){
    for(size_t rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
    if((rnd3 != rnd1) && (rnd3 != rnd2)){

      X[nb1][op_so.id][op_si.id][rnd1][rnd2][rnd3].swap(
      X[nb2][op_so.id][op_si.id][rnd1][rnd2][rnd3]);

    }}}}}// loops random vectors
  }}//loops operators
}
