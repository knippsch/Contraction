/*
 * Perambulator.h
 *
 *  Created on: Apr 13, 2014
 *      Author: werner
 */

#ifndef _PERAMBULATOR_H_
#define _PERAMBULATOR_H_

// TODO: check if they are all necessary. Doesn't matter much though
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "config_utils.h"
#include "EigenVector.h"
#include "global_data.h"
#include "io_utils.h"
#include "propagator_io.h"
#include "quark.h"
#include "RandomVector.h"
#include "typedefs.h"

namespace LapH {

class Perambulator {

public:
  Perambulator ();
  void read_perambulators_from_file (const int config_i);

  // [] operator to directly access the elements of perambulator
  inline const Eigen::MatrixXcd& operator()(const size_t quark_id, 
                                            const size_t rnd_id) const {
    return perambulator[quark_id][rnd_id];
  }

protected:
  array_Xcd_d2_eigen perambulator;

};

} // end of namespace

#endif // _PERAMBULATOR_H__
