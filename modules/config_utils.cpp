/*
 * configs.cpp
 * Source code to config navigation and config mappings
 * Defines: - lookup-tables in 3d (hopping3d)
 *					- mapping from ildg to Eigen Array (map_timeslice_to_eigen)
 *					- building of Laplacian in E4- and C-space (BuildLaplacian)
 *					- gauge transform (transform_ts)
 *					- check for gauge invariance (check_gauge)
 * Created on: Aug 26, 2013
 * Author: christopher helmes
 */

#include "config_utils.h"

static GlobalData * const global_data = GlobalData::Instance();

void hopping3d(int** iup, int** idown){

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();

  int* x0_h = new int[3];
  int* x1_h = new int[3];
  int* x2_h = new int[3];

  int Lt_h = Lx * Ly;

  for ( int x0 = 0; x0 < Lx; ++x0 ) {//loop x0
    x0_h[2] = x0 * Lt_h;
    //negative direction (index at lower boundary)
    if ((x0_h[0] = x0 - 1) < 0) x0_h[0] = Lt_h * (Lx - 1);
    else x0_h[0] *= Lt_h;
    //positive direction (index at upper boundary)
    if ((x0_h[1] = x0 + 1) == Lx) x0_h[1] = 0;
    else x0_h[1] *= Lt_h;

    for ( int x1 = 0; x1 < Lx; ++x1 ) {//loop x1
      x1_h[2] = x1 * Ly;
      //neg. dir.
      if ((x1_h[0] = x1 - 1) < 0) x1_h[0] = Ly * (Lx - 1);
      else x1_h[0] *= Ly;
      //pos. dir.
      if ((x1_h[1] = x1 + 1) == Lx) x1_h[1] = 0;
      else x1_h[1] *= Ly;

      for ( int x2 = 0; x2 < Ly; ++x2 ) {//loop x2
        x2_h[2] = x2;
        //neg. dir.
        if ((x2_h[0] = x2 - 1) < 0) x2_h[0] = Ly -1;
        //pos. dir.
        if ((x2_h[1] = x2 +1) == Ly) x2_h[1] = 0;
        //overall volume index
        int i = x0_h[2] + x1_h[2] + x2_h[2];
        //std::cout << x0 << " " << x1 << " " << x2 << " " << i << std::endl;
        //upwards
        iup[i][0] = x0_h[1] + x1_h[2] + x2_h[2];
        iup[i][1] = x0_h[2] + x1_h[1] + x2_h[2];
        iup[i][2] = x0_h[2] + x1_h[2] + x2_h[1];
        //downwards
        idown[i][0] = x0_h[0] + x1_h[2] + x2_h[2];
        idown[i][1] = x0_h[2] + x1_h[0] + x2_h[2];
        idown[i][2] = x0_h[2] + x1_h[2] + x2_h[0];
      }//end loop x2
    }//end loop x1
  }//end loop x0
  delete x0_h;
  delete x1_h;
  delete x2_h;
}

//Get SU(3)-Matrices from timeslice and sort them into Eigen Array
void map_timeslice_to_eigen(Eigen::Matrix3cd **eigen, double *timeslice) {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();
  const int Vs = Lx * Ly * Lz;
  const int V_TS = global_data->get_V_TS();
  
  //Number of directions
  const int NDIR = 4;
  //Number of colors
  const int NCOL = 3;

  //read in elements
  int el_input = 0;
  for (int z = 0; z < Lz; ++z) {//spatial loops
    for (int y = 0; y < Ly; ++y) {
      for (int x = 0; x < Lx; ++x) {
        for (int mu = 1; mu < 4; ++mu) {//direction loop
          std::complex< double > array[9];
          for (int a = 0; a < 3; ++a) {//colour loops
            for (int b = 0; b < 3; ++b) {
              //timeslice index of real part
              int ind_r = z*V_TS/Lz+y*V_TS/(Lz*Ly)+x*V_TS/(Vs)+
                mu*V_TS/(Vs*NDIR)+a*V_TS/(Lx*Ly*Lz*NDIR*NCOL)
                +b*V_TS/(Vs*NDIR*NCOL*NCOL)+0;
              //timeslice index of imaginary part
              int ind_i = z*V_TS/Lz+y*V_TS/(Lz*Ly)+x*V_TS/(Vs)+
                mu*V_TS/(Vs*NDIR)+a*V_TS/(Vs*NDIR*NCOL)
                +b*V_TS/(Vs*NDIR*NCOL*NCOL)+1;
              std::complex<double> pair(timeslice[ind_r], timeslice[ind_i]);
              //array to be mapped to Eigen Array
              array[3*b+a] = pair;
              ++el_input;
            }
          }
          Eigen::Map<Eigen::Matrix3cd> dummy(array);
          //spatial index
          int ind = z*Ly*Lx+y*Lx+x;
          eigen[ind][mu-1] = dummy;
        }
      }
    }
  }
  //std::cout << el_input << " doubles read in from ildg timeslice " << std::endl;
}

//Displacements

//displacement in one direction i acting to the right
void right_displacement_one_dir(Eigen::Matrix3cd** config, int** iup,
    int** idown, const int dir, Eigen::MatrixXcd& V, Eigen::MatrixXcd& W ) {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();
  const int Vs = Lx * Ly * Lz;

  //Information on Matrix size
  const int num_cols = V.cols();
  const int num_rows = V.rows();

  //Loop over all eigenvectors in 
  for (int i = 0; i < num_cols; ++i ) {
//    std::cout << "eigenvector: " << i << std::endl;
    //storing eigenvector
    Eigen::VectorXcd in(num_rows);
    Eigen::VectorXcd out(num_rows);

    in = V.col(i);

    //Displace eigenvector
    for (int spatial_ind = 0; spatial_ind < Vs; ++spatial_ind) {
      //std::cout << "x: " << spatial_ind << std::endl;
      Eigen::Vector3cd tmp;
      Eigen::Vector3cd quark_up;
      Eigen::Vector3cd quark_down;

      //determine needed indices from lookup tables;
      int up_ind = iup[spatial_ind][dir];
      int down_ind = idown[spatial_ind][dir];

      quark_up = in.segment(3*up_ind,3);
      quark_down = in.segment(3*down_ind,3);
      tmp = 0.5*( (config[spatial_ind][dir] * quark_up) - 
          ((config[down_ind][dir].adjoint()) * quark_down) ); 
      out.segment(3*spatial_ind,3) = tmp;

    }//end spatial loop
    //write displaced eigenvector to W
        W.col(i) = out;
  }//end eigenvector loop

}


