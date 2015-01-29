#ifndef GLOBALDATA_UTILS_H_
#define GLOBALDATA_UTILS_H_

#include <array>
#include <iostream>
#include <string>

#include "global_data_typedefs.h"
#include "typedefs.h"

namespace global_data_utils {

  // functions for input handling
  quark make_quark (const std::string& quark_string);
  void quark_check(quark quarks);
  Operator_list make_operator_list(const std::string& operator_string);

  // functions for lookup tables
  std::array<int, 3> add_p3(const pdg& in1, const pdg& in2);
  int abs_p3(const pdg& in);
  bool compare_quantum_numbers_of_pdg(const pdg& in1, const pdg& in2);
  bool compare_mom_dis_of_pdg(const pdg& in1, const pdg& in2);
  bool compare_index_list(index_IO_1& in1, index_IO_1& in2);
  bool compare_index_list(index_IO_2& in1, index_IO_2& in2);
  bool compare_quantum_numbers_of_pdg(const pdg& in1, const Operators& in2);
  void set_index_corr(vec_pdg_Corr& lookup_corr, vec_pd_VdaggerV& lookup_vdv,
                      vec_pd_rVdaggerVr& lookup_rvdvr);
  void set_index_2pt(const Operators& in1, const Operators& in2,
                     const vec_pdg_Corr& lookup_corr, vec_index_2pt& lookup_2pt);
  void set_index_4pt(const Operators& in1, const Operators& in2, 
                     const Operators& in3, const Operators& in4,
                     const vec_pdg_Corr& lookup_corr, vec_index_4pt& lookup_4pt);
  
} // end of namespace global_data_utils

#endif // GLOBALDATA_UTILS_H_
