/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "global_data_typedefs.h"
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
  int start_config, end_config, delta_config;
  int verbose;
  size_t index_of_unity;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string path_output;
  std::string path_config;
  std::vector<quark> quarks;
  
  std::vector<Operator_list> operator_list;
  Correlator_list correlator_list;

  vec_pdg_Corr lookup_corr;
  vec_index_2pt lookup_2pt;
  vec_index_4pt lookup_4pt;
  vec_index_IO_1 lookup_2pt_IO;
  vec_index_IO_2 lookup_4pt_1_IO;
  vec_index_IO_2 lookup_4pt_2_IO;
  vec_index_IO_1 lookup_4pt_3_IO;
  vec_pd_VdaggerV lookup_vdv;
  vec_pd_rVdaggerVr lookup_rvdvr;
  
  indexlist_1 rnd_vec_1pt;
  indexlist_2 rnd_vec_2pt;
  indexlist_3 rnd_vec_3pt;
  indexlist_4 rnd_vec_4pt;

  void init_lookup_tables();

  void input_handling(const std::vector<std::string>& quark_configs,
                      const std::vector<std::string>& operator_list_configs,
                      const std::vector<std::string>& correlator_list_configs);

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
  inline std::vector<Operator_list>& get_operator_list() {
    return operator_list;
  }
  inline Correlator_list& get_correlator_list() {
    return correlator_list;
  }
  inline const vec_pdg_Corr& get_lookup_corr() {
    return lookup_corr;
  }
  inline const vec_index_2pt& get_lookup_2pt_trace() {
    return lookup_2pt;
  }
  inline const vec_index_4pt& get_lookup_4pt_trace() {
    return lookup_4pt;
  }
  inline const vec_index_IO_1& get_lookup_2pt_IO() {
    return lookup_2pt_IO;
  }
  inline const vec_index_IO_2& get_lookup_4pt_1_IO() {
    return lookup_4pt_1_IO;
  }
  inline const vec_index_IO_2& get_lookup_4pt_2_IO() {
    return lookup_4pt_2_IO;
  }
  inline const vec_index_IO_1& get_lookup_4pt_3_IO() {
    return lookup_4pt_3_IO;
  }
  inline const indexlist_1& get_rnd_vec_1pt() {
    return rnd_vec_1pt;
  }
  inline const indexlist_2& get_rnd_vec_2pt() {
    return rnd_vec_2pt;
  }
  inline const indexlist_3& get_rnd_vec_3pt() {
    return rnd_vec_3pt;
  }
  inline const indexlist_4& get_rnd_vec_4pt() {
    return rnd_vec_4pt;
  }
  inline const size_t get_index_of_unity() {
    return index_of_unity;
  }
  inline const vec_pd_VdaggerV get_lookup_VdaggerV() {
    return lookup_vdv;
  }
  inline const vec_pd_rVdaggerVr get_lookup_rVdaggerVr() {
    return lookup_rvdvr;
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

#endif /* GLOBALDATA_H_ */
