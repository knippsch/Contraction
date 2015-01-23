#if 0
/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

//#ifndef GLOBALDATA_H_
//#define GLOBALDATA_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "quark.h"
#include "typedefs.h"

class GlobalData {

private:
  //! A pointer on the class itself
  static GlobalData* instance_;
  //! globally accessible data
  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_rnd_vec;
  int number_of_inversions;
  int number_of_max_mom;
  int max_mom_in_one_dir;
  int number_of_momenta;
  int dirac_min;
  int dirac_max;
  int displ_min;
  int displ_max;
  size_t number_of_displ;
  size_t number_of_dirac;
  int start_config, end_config, delta_config;
  int verbose;
  size_t number_of_operators;
  size_t number_of_displ_gamma;
  size_t number_of_momentum_squared;
  size_t index_of_unity;
  size_t number_of_VdaggerV;
  size_t number_of_rVdaggerVr;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string path_output;
  std::string path_config;
  std::vector<quark> quarks;
  std::vector<int> momentum_squared;
  void quark_input_data_handling (const std::vector<std::string> quark_configs);

////////////////////////////////////////////////////////////////////////////////
  //TODO: Clean that up

  //dirac structure hardcoded
  std::vector<size_t> dg;

  vec_pdg_Corr op_Corr;
  
  vec_pd_VdaggerV op_VdaggerV;
  vec_pd_rVdaggerVr op_rVdaggerVr;

  vec_pdg_C2 op_C2;
  vec_pdg_C4 op_C4;

  indexlist_1 rnd_vec_C1;
  indexlist_2 rnd_vec_C2;
  indexlist_3 rnd_vec_C3;
  indexlist_4 rnd_vec_C4;

  void init_from_infile();
  void set_Corr();
  void set_C2();
  void set_C4();
  void set_rnd_vec_C1();
  void set_rnd_vec_C2();
  void set_rnd_vec_C3();
  void set_rnd_vec_C4();
////////////////////////////////////////////////////////////////////////////////


public:

  static GlobalData* Instance ();

  void read_parameters(int ac, char* av[]);

  inline std::string get_name_lattice() {
    return name_lattice;
  }
  inline std::string get_output_path() {
    return path_output;
  }
  inline std::string get_config_path() {
    return path_config;
  }
  inline int get_Lx () {
    return Lx;
  }
  inline int get_Ly () {
    return Ly;
  }
  inline int get_Lz () {
    return Lz;
  }
  inline int get_Lt () {
    return Lt;
  }
  inline int get_dim_row () {
    return dim_row;
  }
  inline int get_V_TS () {
    return V_TS;
  }
  inline int get_V_for_lime () {
    return V_for_lime;
  }
  inline int get_number_of_inversions () {
    return number_of_inversions;
  }
  inline int get_number_of_rnd_vec () {
    return number_of_rnd_vec;
  }
  inline int get_number_of_max_mom () {
    return number_of_max_mom;
  }
  inline int get_max_mom_in_one_dir () {
    return max_mom_in_one_dir;
  }
  inline int get_number_of_momenta() {
    return momentum_squared.size();
  }
  inline const std::vector<size_t>& get_displ_gamma() {
    return dg;
  }
  inline int get_dirac_min() {
    return dirac_min;
  }
  inline int get_dirac_max() {
    return dirac_max;
  }
  inline int get_displ_min() {
    return displ_min;
  }
  inline int get_displ_max() {
    return displ_max;
  }
  inline int get_start_config () {
    return start_config;
  }
  inline int get_end_config () {
    return end_config;
  }
  inline int get_delta_config () {
    return delta_config;
  }
  inline int get_number_of_eigen_vec() {
    return number_of_eigen_vec;
  }
  inline int get_verbose() {
    return verbose;
  }
  inline size_t get_number_of_operators() {
    return number_of_operators;
  }
  inline size_t get_number_of_displ_gamma() {
    return number_of_displ_gamma;
  }
  inline size_t get_number_of_displ() {
    return number_of_displ;
  }
  inline size_t get_number_of_dirac() {
    return number_of_dirac;
  }
  inline size_t get_number_of_momentum_squared() {
    return number_of_momentum_squared;
  }
  inline std::string get_path_eigenvectors() {
    return path_eigenvectors;
  }
  inline std::string get_name_eigenvectors() {
    return name_eigenvectors;
  }
  inline std::string get_path_perambulators() {
    return path_perambulators;
  }
  inline std::string get_name_perambulators() {
    return name_perambulators;
  }
  inline std::vector<quark> get_quarks() {
    return quarks;
  }
  inline std::vector<int> get_momentum_squared() {
    return momentum_squared;
  }
  inline const vec_pdg_Corr& get_op_Corr() {
    return op_Corr;
  }
  inline const vec_pdg_C2& get_op_C2() {
    return op_C2;
  }
  inline const vec_pdg_C4& get_op_C4() {
    return op_C4;
  }
  inline const indexlist_1& get_rnd_vec_C1() {
    return rnd_vec_C1;
  }
  inline const indexlist_2& get_rnd_vec_C2() {
    return rnd_vec_C2;
  }
  inline const indexlist_3& get_rnd_vec_C3() {
    return rnd_vec_C3;
  }
  inline const indexlist_4& get_rnd_vec_C4() {
    return rnd_vec_C4;
  }
  inline const size_t get_number_of_VdaggerV() {
    return number_of_VdaggerV;
  }
  inline const size_t get_number_of_rVdaggerVr() {
    return number_of_rVdaggerVr;
  }
  inline const size_t get_index_of_unity() {
    return index_of_unity;
  }
  inline const vec_pd_VdaggerV get_op_VdaggerV() {
    return op_VdaggerV;
  }
  inline const vec_pd_rVdaggerVr get_op_rVdaggerVr() {
    return op_rVdaggerVr;
  }

  //! All con/de-structors are protected to assure that only one instance exists
  //! at once. DO NOT CHANGE!!
protected:
  GlobalData () {
  }
  GlobalData (const GlobalData& other) {
  }
  virtual ~GlobalData () {
  }

};

//#endif /* GLOBALDATA_H_ */
#endif
