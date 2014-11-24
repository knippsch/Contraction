#include "CrossOperator.h"  

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

LapH::CrossOperator::CrossOperator(const size_t number) : X(number) {

  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // TODO: must be changed in GlobalData {
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

  for(auto& xx : X){
    xx.resize(boost::extents[nb_mom][nb_mom][nb_dir][nb_dir]
                            [nb_rnd][nb_rnd][nb_rnd]);

    std::fill(xx.data(), xx.data() + xx.num_elements(), 
                Eigen::MatrixXcd(4 * dilE, 4 * dilE));
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::construct(const BasicOperator& basic, 
                                    const VdaggerV& vdaggerv, const size_t nb,
                                    const int t_source, const int t_sink){

  const int Lt = global_data->get_Lt();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t dilT = quarks[0].number_of_dilution_T;
  // TODO: must be changed in GlobalData {
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }

  const std::array<double, 4> bla = {{1., 1., -1., -1.}};
  size_t tu, td, t2;
  if(nb == 0){
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

  for(size_t p_so = 0; p_so < nb_mom; p_so++){
  for(size_t p_si = 0; p_si < nb_mom; p_si++){
    for(size_t d_so = 0; d_so < nb_dir; d_so++){
    for(size_t d_si = 0; d_si < nb_dir; d_si++){
      #pragma omp parallel for collapse(2) schedule(dynamic)
      for(size_t rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
      for(size_t rnd2 = 0; rnd2 < nb_rnd; ++rnd2){
      if(rnd2 != rnd1){
      for(size_t rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
      if((rnd3 != rnd1) && (rnd3 != rnd2)){
        for(size_t col = 0; col < 4; col++){
        for(size_t row = 0; row < 4; row++){

            X[nb][p_so][p_si][d_so][d_si][rnd1][rnd2][rnd3]
                                    .block(row*dilE, col*dilE, dilE, dilE) = 

            bla[col] * //TODO: Dirac by hand is not good! 

            basic.get_operator(t_source, tu, td, d_so, p_so, rnd1, rnd2)
                      .block(row*dilE, col*dilE, dilE, dilE) *

            vdaggerv.return_rvdaggervr(p_si, t2, d_si, rnd2, rnd3)
                      .block(0, col*dilE, dilE, dilE);

        }}// loops over col and row end here
      }}}}}// loops random vectors
    }}// loops dirac indices
  }}// loops momenta
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::swap(const size_t nb1, const size_t nb2){
  
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  // TODO: must be changed in GlobalData {
  std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();
  // TODO: }
  // TODO: Think about for each loop
  #pragma omp parallel for collapse(6) schedule(dynamic)
  for(size_t p_so = 0; p_so < nb_mom; p_so++){
  for(size_t p_si = 0; p_si < nb_mom; p_si++){
    for(size_t d_so = 0; d_so < nb_dir; d_so++){
    for(size_t d_si = 0; d_si < nb_dir; d_si++){
      for(size_t rnd1 = 0; rnd1 < nb_rnd; ++rnd1){
      for(size_t rnd2 = 0; rnd2 < nb_rnd; ++rnd2){
      if(rnd2 != rnd1){
      for(size_t rnd3 = 0; rnd3 < nb_rnd; ++rnd3){
      if((rnd3 != rnd1) && (rnd3 != rnd2)){

        X[nb1][p_so][p_si][d_so][d_si][rnd1][rnd2][rnd3].swap(
        X[nb2][p_so][p_si][d_so][d_si][rnd1][rnd2][rnd3]);

      }}}}}// loops random vectors
    }}// loops dirac indices
  }}// loops momenta
}
