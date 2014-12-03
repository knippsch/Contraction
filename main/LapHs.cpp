//============================================================================
// Name        : LapHs.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : stochastic LapH code
//============================================================================

#include <iostream>

#include "Correlators.h"
#include "GlobalData.h"

#include "omp.h"

int main (int ac, char* av[]) {

  // initialization of OMP paralization
  // TODO: number of thread must be set via inputfile
  Eigen::initParallel();
  omp_set_num_threads(2);
  Eigen::setNbThreads(1); // parallel eigen makes everything slower

  // reading in global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  LapH::Correlators corr;

  // ***************************************************************************
  // Loop over all configurations **********************************************
  // ***************************************************************************
  const int end_config = global_data->get_end_config();
  const int delta_config = global_data->get_delta_config();
  const int start_config = global_data->get_start_config();

  for(size_t config_i = start_config; config_i <= end_config; config_i +=
      delta_config){
    std::cout << "\nprocessing configuration: " << config_i << "\n\n";
    corr.compute_correlators(config_i);
  }

}

