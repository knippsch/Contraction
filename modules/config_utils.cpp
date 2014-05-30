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
          Eigen::Map<Eigen::Matrix3cd> dummy (array);
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


//HYP-Smearing and helper functions


//storage map for dec_timeslice
static int decor (int dir, int smear_plane) {
  int ret;
  if (dir == 0) {
    if (smear_plane == 1) ret = 0;
    else ret = 3; 
  }
  else if (dir == 1) {
    if (smear_plane == 0) ret = 1;
    else ret = 4;
  }
  else {
    if (smear_plane == 0) ret = 2;
    else ret = 5;
  }
  return ret;
}


static std::array< int, 2 > get_dirs(int mu) {
  std::array<int,2> dirs;  
  if (mu == 0) {
    dirs.at(0) = 1;
    dirs.at(1) = 2;
  }
  else if (mu == 1) {
    dirs.at(0) = 0;
    dirs.at(1) = 2; 
  }
  else if (mu == 2) {
    dirs.at(0) = 0;
    dirs.at(1) = 1; 
  }
  
  return(dirs);
}

//project a 3x3matrix back to SU3 form
static Eigen::Matrix3cd proj_to_su3_imp(Eigen::Matrix3cd& in){
  //avoid possible aliasing issues:
  Eigen::Matrix3cd out = ( (in.adjoint() * in).sqrt() ).inverse();
  in = in*out;
  std::complex<double> det = 1./pow(in.determinant(),1./3.);
  in *= det;
  return in;
}

void smearing_hyp(int** iup, int** idown, Eigen::Matrix3cd **eigen_timeslice, 
    double alpha_1, double alpha_2, int iter) {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();
  const int Vs = Lx * Ly * Lz;

  // temporal timeslice twice the size for decorated links 
  Eigen::Matrix3cd **dec_timeslice = new Eigen::Matrix3cd *[Vs];
  for (auto vol = 0; vol < Vs; ++vol) {
    dec_timeslice[vol] = new Eigen::Matrix3cd[6];
  }

  // temporal timeslice from decorated links
  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[Vs]; 
  for ( auto i = 0; i < Vs; ++i ) {
    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  }
  
  
  // temporal integers holding directions for decorated smearing
  int mu, nu, eta;
  for (int run = 0; run < iter; ++run) {
    // calculate inner staple from original timeslice, store in dec_timeslice 
    // each link can get smeared in two planes
    for (auto vol = 0; vol < Vs; ++vol) {
      for (auto dir = 0; dir < 3; ++dir) {
        //inner staple
        Eigen::Matrix3cd inner_staple = Eigen::Matrix3cd::Zero();
        std::array< Eigen::Matrix3cd, 2 > tmp_staples;
        std::array<int, 2> perpendics = get_dirs( dir );
        for (auto it_perp_dir = perpendics.begin(); 
            it_perp_dir != perpendics.end(); ++it_perp_dir ) {
          int perp_dir = *it_perp_dir;
          //up-type smearing_indices
          mu = iup[vol][perp_dir];
          nu = iup[mu][dir];
          eta = iup[vol][dir];

          //product of up matrices
          inner_staple = eigen_timeslice[vol][perp_dir] * 
              ( eigen_timeslice[mu][dir] * 
              ( eigen_timeslice[eta][perp_dir].adjoint() ) ); 
          //down-type smearing indices
          mu = idown[vol][perp_dir];
          nu = iup[mu][dir];

          //eta is same endpoint no adjoint necessary here
          //product of down matrices
          inner_staple += ( eigen_timeslice[mu][perp_dir].adjoint() ) *
              ( eigen_timeslice[mu][dir] * eigen_timeslice[nu][perp_dir] );

          //Careful placement of decorated links in dec_timeslices:
          //dir=0 has placement in dec_dir = 0 (smeared in 1 plane)
          //                       dec_dir = 3 (smeared in 2 plane)
          //dir=1 has placement in dec_dir = 1 (smeared in 0 plane)
          //                       dec_dir = 4 (smeared in 2 plane)
          //dir=2 has placement in dec_dir = 2 (smeared in 0 plane)
          //                       dec_dir = 5 (smeared in 1 plane)
          Eigen::Matrix3cd stac = ( eigen_timeslice[vol][dir] *
              (1-alpha_2) ) +  ( inner_staple * alpha_2/2.);  
          int n_el = it_perp_dir - perpendics.begin();
          tmp_staples.at(n_el) = proj_to_su3_imp(stac);
          //tmp_staples.at(n_el) = stac;//without SU(3) projection
        }

        //staple link in direction dir in non participating and negative directions
        dec_timeslice[vol][dir] = tmp_staples.at(0);
        dec_timeslice[vol][dir+3] = tmp_staples.at(1);
      }
    }
    //just debugging
    //if(run == 0) std::cout << "Fat link a x,y,z (11,19,29): \n" << dec_timeslice[11*L2*L1+19*L1+29][1]<< "\n\n";
    //calculate outer staple from dec_timeslice as modified ape-smearing

      for (int i = 0; i < Vs; ++i) {
        for (int dir = 0; dir < 3; ++dir) {
          Eigen::Matrix3cd outer_staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link

          //filling each element of smearer using 3 links
          //debugging link
          Eigen::Matrix3cd outer_staple_test;
          for (int not_dir = 0; not_dir < 3; ++not_dir) {
            if (dir != not_dir) {

              //calculate plane in which was smeared
              int plane = ( (dir+1) ^ (not_dir+1) ) - 1;
              mu = iup[i][not_dir];
              nu = iup[mu][dir];
              eta = idown[nu][not_dir];

              //Staples in positive direction
              //replace directions by appropriate decor 
              int a,b;
              a = decor(not_dir,plane);
              b = decor(dir,plane);
              outer_staple += dec_timeslice[i][a]*
                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint()));
/*           if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" << dec_timeslice[i][a]*
                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint())) << "\n\n";
                */
              mu = idown[i][not_dir];
              nu = iup[mu][dir];
            
              //Staples in negative direction
              outer_staple += (dec_timeslice[mu][a].adjoint())*
                (dec_timeslice[mu][b]*dec_timeslice[nu][a]); //structure has to be a_dag, b, a
  /*         if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" <<(dec_timeslice[mu][a].adjoint())*
                (dec_timeslice[mu][b]*dec_timeslice[nu][a]) << "\n\n";
                */
 
            }
          }
          eigen_timeslice_ts[i][dir] = (eigen_timeslice[i][dir] * (1.-alpha_1)) + (outer_staple * alpha_1/4.);
        }
      }
      for ( auto i = 0; i < Vs; ++i ) {
        for ( auto mu = 0; mu < 3; ++mu) {
          eigen_timeslice[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
          //eigen_timeslice[i][mu] = eigen_timeslice_ts[i][mu];//without SU(3)-projection
        }
      }
  }
  //clean up
  for (int k = 0; k < 3; ++k) {
    delete[] dec_timeslice[k];
    delete[] dec_timeslice[k+3];
    delete[] eigen_timeslice_ts[k];
  } 
  delete dec_timeslice;
  delete eigen_timeslice_ts;
  //printf("Timeslice successfully HYP-smeared with (a_1, a_2, iterations): %f, %f, %d \n",
      //alpha_1, alpha_2, iter);
}

