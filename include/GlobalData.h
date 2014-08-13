/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <algorithm>
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

#include "quark.h"

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
  int start_config, end_config, delta_config;
  int verbose;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string path_output;
  std::vector<quark> quarks;
  std::vector<int> momentum_squared;
  void quark_input_data_handling (const std::vector<std::string> quark_configs);


public:
  static GlobalData* Instance ();

  void read_parameters(int ac, char* av[]);

  inline std::string get_name_lattice() {
    return name_lattice;
  }
  inline std::string get_output_path() {
    return path_output;
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
