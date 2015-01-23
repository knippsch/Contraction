#include "global_data.h"
#include "global_data_utils.h"

namespace {

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
static void copy_quantum_numbers(const pdg& in, std::array<int, 6>& out){

  out[0] = in.dis3[0];
  out[1] = in.dis3[1];
  out[2] = in.dis3[2];
  out[3] = in.p3[0];
  out[4] = in.p3[1]; 
  out[5] = in.p3[2];
}


} // end of unnamed namespace

namespace global_data_utils {

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// function that compares two pdg structs and checks if the corresponding 
// entries of lookup_corr coincide
bool compare_quantum_numbers_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3) && 
      (in1.gamma == in2.gamma))
    return true;
  else
    return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
bool compare_mom_dis_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3))
    return true;
  else
    return false;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// function that compares a pdg struct with an operator as defined by the input
// file and checkes if the quantum numbers of pdg are contained in the physical
// situation described by the operator
bool compare_quantum_numbers_of_pdg(const pdg& in1, const Operators& in2){

  if( in1.gamma == in2.gammas &&
      in1.dis3 == in2.dil_vec){
    for(const auto& mom_vec : in2.mom_vec){
      for(auto mom : mom_vec){
        if(in1.p3 == mom){
          return true;
          //TODO: is that safer or just spam?
          break;
        }
      }
    }
  }

  return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
bool compare_index_list(index_IO_1& in1, index_IO_1& in2) {

  if(in1.index_pt.size() != in2.index_pt.size())
    return false;

  auto it1 = in1.index_pt.begin();
  auto it2 = in2.index_pt.begin();
  while(it1 != in1.index_pt.end()){
    if(*it1 != *it2){
      return false;
      break;
    }
  }

  return true;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
bool compare_index_list(index_IO_2& in1, index_IO_2& in2) {

  if(in1.index_pt.size() != in2.index_pt.size())
    return false;

  auto it1 = in1.index_pt.begin();
  auto it2 = in2.index_pt.begin();
  while(it1 != in1.index_pt.end()){
    if(*it1 != *it2){
      return false;
      break;
    }
  }

  return true;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void set_index_corr(vec_pdg_Corr& lookup_corr, vec_pd_VdaggerV& lookup_vdv,
                    vec_pd_rVdaggerVr& lookup_rvdvr) {

  // first number is the operator id, the next three is the displacement vector
  // and the last three the momentum vector
  std::vector<std::array<int, 6> > rvdaggervr_qu_nb;
  std::vector<std::array<int, 6> > vdaggerv_qu_nb;
  size_t counter_rvdvr = 0;
  size_t counter_vdv = 0;
  for(auto& op : lookup_corr){
    std::array<int, 6> write;
    if(op.id != 0){
      copy_quantum_numbers(op, write);
      // ######################################################################
      // check if quantum numbers are already stored in rvdaggervr_qu_nb
      bool is_known_rvdvr = false;
      size_t fast_counter_rvdvr = 0;// this gives the Op id if QN are duplicate
      for(const auto& rvdvr : rvdaggervr_qu_nb){
        if(rvdvr == write){
          is_known_rvdvr = true;
          break;
        }
        fast_counter_rvdvr++;
      }
      if(!is_known_rvdvr){ // setting the unknown quantum numbers
        op.id_rvdvr = counter_rvdvr;
        counter_rvdvr++;
        rvdaggervr_qu_nb.push_back(write);
      }
      else
        op.id_rvdvr = fast_counter_rvdvr;

      // ######################################################################
      // check if quantum numbers are already stored in vdaggerv_qu_nb
      bool is_known_vdv = false;
      op.first_vdv = false;
      op.negative_momentum = false;
      size_t fast_counter_vdv = 0;// this gives the Op id if QN are duplicate
      // first check for duplicate quantum numbers
      for(const auto& vdv : vdaggerv_qu_nb){
        if(vdv == write){
          is_known_vdv = true;
          break;
        }
        fast_counter_vdv++;
      }
      if(!is_known_vdv){ // second check for complex conjugate momenta
        fast_counter_vdv = 0;
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;
        for(const auto& vdv : vdaggerv_qu_nb){
          if(vdv == write){
            is_known_vdv = true;
            break;
          }
          fast_counter_vdv++;
        }
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;

        // case quantum numbers are unknown -> create new entry for 
        // vdaggerv_qu_nb
        if(!is_known_vdv){
          op.id_vdv = counter_vdv;
          vdaggerv_qu_nb.push_back(write);
          counter_vdv++;
          op.first_vdv = true;
        }
        // case quantum numbers are unknown, but opposite momente exist. 
        // VdaggerV can be obtained from the negative momentum.
        else{
          op.negative_momentum = true;
          op.id_vdv = fast_counter_vdv;
        }
      }
      // case same quantum numbers already exist. 
      // TODO: I don't think that works if several negative momenta exist.
      else{
        op.negative_momentum = lookup_corr[fast_counter_vdv].negative_momentum;
        op.id_vdv = fast_counter_vdv;
      }
    }
    else{ // setting the very first entry
      copy_quantum_numbers(op, write);
      rvdaggervr_qu_nb.push_back(write);
      vdaggerv_qu_nb.push_back(write);
      op.id_vdv = counter_vdv;
      op.id_rvdvr = counter_rvdvr;
      op.first_vdv = true;
      counter_vdv++;
      counter_rvdvr++;
    }
  }

  // setting the lookuptables to be able to reconstruct the quantum numbers
  // when computing VdaggerV and rVdaggerVr
  lookup_vdv.resize(vdaggerv_qu_nb.size());
  lookup_rvdvr.resize(rvdaggervr_qu_nb.size());

  size_t index = 0;
  for(auto& op_vdv : lookup_vdv){
    op_vdv.id = index;
    for(const auto& op : lookup_corr){
      if(index == op.id_vdv)
        op_vdv.index = op.id;
    }
    index++;
  }

  index = 0;
  for(auto& op_rvdvr : lookup_rvdvr){
    op_rvdvr.id = index;
    for(const auto& op : lookup_corr){
      if(index == op.id_vdv){
        op_rvdvr.index = op.id;
        // if the momentum was only build for the negative in VdaggerV, the 
        // adjoint has been taken
        if(op.negative_momentum == true)
          op_rvdvr.adjoint = true;
        else
          op_rvdvr.adjoint = false;
      }
    }
    index++;
  }

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void set_index_2pt(const Operators& in1, const Operators& in2, 
                   const vec_pdg_Corr& lookup_corr, vec_index_2pt& lookup_2pt) {

  index_2pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      write.index_Q2 = op1.id;
      write.index_Corr = op2.id;

      lookup_2pt.push_back(write);
    }} //loops over sink end here
  }} //loops over source end here

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void set_index_4pt(const Operators& in1, const Operators& in2, 
                   const Operators& in3, const Operators& in4,
                   const vec_pdg_Corr& lookup_corr, vec_index_4pt& lookup_4pt) {

  index_4pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      for(const auto& op3 : lookup_corr){
      if(compare_quantum_numbers_of_pdg(op3, in3)){
        for(const auto& op4 : lookup_corr){
        if(compare_quantum_numbers_of_pdg(op4, in4)){
          write.index_Q2[0]   = op1.id;
          write.index_Corr[0] = op2.id;
          write.index_Q2[1]   = op3.id;
          write.index_Corr[1] = op4.id;

          lookup_4pt.push_back(write);
        }}
      }} //loops over sink end here
    }}
  }} //loops over source end here

}


} // end of namespace global_data_utils
