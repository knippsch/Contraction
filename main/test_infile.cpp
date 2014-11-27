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
#include "OperatorStructure.h"

int main (int ac, char* av[]) {

  // reading in global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  const int start_config = global_data->get_start_config();

  std::cout << "\nprocessing configuration: " << start_config << "\n\n";

  std::vector<LapH::pdg> op_Corr;
  std::vector<LapH::pdg_C2> op_C2;

  LapH::init_from_infile(op_Corr, op_C2);


//  for(auto& i : op_Corr) std::cout << "id " << i.id << std::endl;
//  std::cout << "p " << std::endl;
//  for(auto& i : op_Corr) std::cout << i.p[0] << i.p[1] << i.p[2] << std::endl;
//  std::cout << "dis " << std::endl;
//  for(auto& i : op_Corr) std::cout << i.dis[0] << i.dis[1] << i.dis[2] << std::endl;
//  std::cout << "gamma " << std::endl;
//  for(auto& i : op_Corr) std::cout << i.gamma[0] << i.gamma[1] << i.gamma[2] << i.gamma[3] << std::endl;

  std::cout << "p " << std::endl;
  for(auto& i : op_C2) for(auto& mom : i.index) std::cout << i.p_sq << i.dg_so 
    << i.dg_si << "  " << mom[0] << "  " << mom[1] << std::endl;




  std::cout << op_C2.size() << std::endl;
  std::cout << op_Corr.size() << std::endl;
  
  return 0;
}

