#ifndef CORRELATORS_H_
#define CORRELATORS_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "GlobalData.h"
#include "BasicOperator.h"
#include "Perambulator.h"
#include "typedefs.h"

#include "omp.h"

namespace LapH {

class Correlators{

public:
  Correlators();
  ~Correlators() {};

  void build_everything(const size_t config_i);

private:
   BasicOperator basic;
   array_cd_d5 C4_mes;
   array_cd_d6 C2_mes;
   array_cd_d10 Corr;

   void compute_meson_corr(const int t_source, const int t_sink);
   void compute_meson_4pt_cross(array_Xcd_d7_eigen& X, array_Xcd_d7_eigen& Y, 
                                const int t_source, const int t_sink);

};

}// end of namespace

#endif // CORRELATORS_H_
