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
  const std::vector<quark> quarks = global_data->get_quarks();
  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_C3();
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
#pragma omp single
{     
  for(const auto& op : op_C4){
  for(const auto& i : op.index){
    size_t id_so = i[2+nb];
    size_t id_si = i[1-nb];

    #pragma omp task shared (op, i)
    for(auto& rnd_it : rnd_vec_index) {
      compute_X(basic, id_si, 
                basic.get_operator(t_source, tu, td, id_so, rnd_it[0], rnd_it[1]),
                vdaggerv.return_rvdaggervr(id_si, t2, rnd_it[1], rnd_it[2]),
                X[nb][id_so][id_si][rnd_it[0]][rnd_it[1]][rnd_it[2]]);

//      compute_X(basic, vdaggerv, id_so, id_si, rnd1, rnd2, rnd3, t_source, tu, 
//          td, t2, X[nb]);

    } // loop over random vectors
  }}// loops operators
}

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::compute_X(const BasicOperator& basic, 
                                    const size_t id_si, 
                                    const Eigen::MatrixXcd& Q2, 
                                    const Eigen::MatrixXcd& VdaggerV,
                                    Eigen::MatrixXcd& X) {

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;

  for(size_t block = 0; block < 4; block++){

    cmplx value = 1;
    basic.value_dirac(id_si, block, value);

    X.block(0, block*dilE, 4*dilE, dilE) = 
      value * Q2.block(0, block*dilE, 4*dilE, dilE) * 
      VdaggerV.block(0, basic.order_dirac(id_si, block)*dilE, dilE, dilE);

    }// loop block ends here

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::swap(const size_t nb1, const size_t nb2){
  
  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_C3();
  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();
  // TODO: }
  // TODO: Think about for each loop

  // omp parallel for cannot handle autoloops,
  // therefore this workaround is implemented
  #pragma omp parallel
  #pragma omp single
  {
    for(auto& op_so : op_Corr){
    for(auto& op_si : op_Corr){
      #pragma omp task shared(op_so, op_si)
      for(auto& rnd_it : rnd_vec_index) {
        X[nb1][op_so.id][op_si.id][rnd_it[0]][rnd_it[1]][rnd_it[2]].swap(
        X[nb2][op_so.id][op_si.id][rnd_it[0]][rnd_it[1]][rnd_it[2]]);

      } // random vector loop
    }}//loops operators
  }
}
