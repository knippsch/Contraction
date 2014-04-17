/*
 * BasicOperator.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#include "BasicOperator.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();


namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);

static void create_gamma (
    Eigen::SparseMatrix<std::complex<double> >* const gamma, const int i) {

  try{
    switch(i) {
    case 0: // gamm_0 (time component)
      gamma[0].insert(0, 2) = -1.;
      gamma[0].insert(1, 3) = -1.;
      gamma[0].insert(2, 0) = -1.;
      gamma[0].insert(3, 1) = -1.;
      break;
    case 1: // gamma_1
      gamma[1].insert(0, 3) = -I;
      gamma[1].insert(1, 2) = -I;
      gamma[1].insert(2, 1) = I;
      gamma[1].insert(3, 0) = I;
      break;
    case 2: // gamma_2
      gamma[2].insert(0, 3) = -1.;
      gamma[2].insert(1, 2) = 1.;
      gamma[2].insert(2, 1) = 1.;
      gamma[2].insert(3, 0) = -1.;
      break;
    case 3: // gamma_3
      gamma[3].insert(0, 2) = -I;
      gamma[3].insert(1, 3) = I;
      gamma[3].insert(2, 0) = I;
      gamma[3].insert(3, 1) = -I;
      break;
    case 4: // unity
      gamma[4].insert(0, 0) = 1.;
      gamma[4].insert(1, 1) = 1.;
      gamma[4].insert(2, 2) = 1.;
      gamma[4].insert(3, 3) = 1.;
      break;
    case 5: // gamma_5
      gamma[5].insert(0, 0) = 1.;
      gamma[5].insert(1, 1) = 1.;
      gamma[5].insert(2, 2) = -1.;
      gamma[5].insert(3, 3) = -1.;
      break;
    case 6: // gamma_0 gamma_5
      gamma[6].insert(0, 2) = 1.;
      gamma[6].insert(1, 3) = 1.;
      gamma[6].insert(2, 0) = -1.;
      gamma[6].insert(3, 1) = -1.;
      break;
    case 7: // gamma_1 gamma_5
      gamma[7].insert(0, 3) = I;
      gamma[7].insert(1, 2) = I;
      gamma[7].insert(2, 1) = I;
      gamma[7].insert(3, 0) = I;
      break;
    case 8: // gamma_2 gamma_5
      gamma[8].insert(3, 0) = 1.;
      gamma[8].insert(2, 1) = -1.;
      gamma[8].insert(1, 2) = 1.;
      gamma[8].insert(0, 3) = -1.;
      break;
    case 9: // gamma_3 gamma_5
      gamma[9].insert(0, 2) = I;
      gamma[9].insert(1, 3) = -I;
      gamma[9].insert(2, 0) = I;
      gamma[9].insert(3, 1) = -I;
      break;
    case 10: // gamma_0 gamma_1
      gamma[10].insert(0, 1) = -I;
      gamma[10].insert(1, 0) = -I;
      gamma[10].insert(2, 3) = I;
      gamma[10].insert(3, 2) = I;
      break;
    case 11: // gamma_0 gamma_2
      gamma[11].insert(0, 1) = -1.;
      gamma[11].insert(1, 0) = 1.;
      gamma[11].insert(2, 3) = 1.;
      gamma[11].insert(3, 2) = -1.;
      break;
    case 12: // gamma_0 gamma_3
      gamma[12].insert(0, 0) = -I;
      gamma[12].insert(1, 1) = I;
      gamma[12].insert(2, 2) = I;
      gamma[12].insert(3, 3) = -I;
      break;
    case 13: // gamma_1 gamma_2
      gamma[13].insert(0, 0) = I;
      gamma[13].insert(1, 1) = -I;
      gamma[13].insert(2, 2) = I;
      gamma[13].insert(3, 3) = -I;
      break;
    case 14: // gamma_1 gamma_3
      gamma[14].insert(0, 1) = -1.;
      gamma[14].insert(1, 0) = 1.;
      gamma[14].insert(2, 3) = -1.;
      gamma[14].insert(3, 2) = 1.;
      break;
    case 15: // gamma_2 gamma_3
      gamma[15].insert(0, 1) = I;
      gamma[15].insert(1, 0) = I;
      gamma[15].insert(2, 3) = I;
      gamma[15].insert(3, 2) = I;
      break;
    case 16: // gamma_2 gamma_0 gamma_5
      gamma[16].insert(0, 1) = 1.;
      gamma[16].insert(1, 0) = -1.;
      gamma[16].insert(2, 3) = 1.;
      gamma[16].insert(3, 2) = -1.;
      break;
    case 17: // nothing yet
      gamma[17].insert(0, 0) = 1.;
      gamma[17].insert(1, 1) = 1.;
      gamma[17].insert(2, 2) = 1.;
      gamma[17].insert(3, 3) = 1.;
      break;
    case 18: // nothing yet
      gamma[18].insert(0, 0) = 1.;
      gamma[18].insert(1, 1) = 1.;
      gamma[18].insert(2, 2) = 1.;
      gamma[18].insert(3, 3) = 1.;
      break;
    case 19: // nothing yet
      gamma[19].insert(0, 0) = 1.;
      gamma[19].insert(1, 1) = 1.;
      gamma[19].insert(2, 2) = 1.;
      gamma[19].insert(3, 3) = 1.;
      break;
    default:
      printf("Dirac component %d not found\n", i);
      exit(0);
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static void create_momenta (std::complex<double>** momentum) {

  try{
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();
    // helper variables for momenta
    const double px = 2. * M_PI / (double) Lx;
    const double py = 2. * M_PI / (double) Ly;
    const double pz = 2. * M_PI / (double) Lz;
    int p = 0;
    // running over all momentum components
    for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
      for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
        for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
          // running over all lattice points
          for(int x = 0; x < Lx; ++x){
            const int xH = x * Ly * Lz; // helper variable
            const int ipxH = ipx * px * x; // helper variable
            for(int y = 0; y < Ly; ++y){
              const int xHyH = xH + y * Lz; // helper variable
              const int ipxHipzH = ipxH + ipy * py * y; // helper variable
              for(int z = 0; z < Lz; ++z){
                //ipz=1;
                momentum[p][xHyH + z] = exp(-I * (ipxHipzH + ipz * pz * z));
                //std::cout << "mom = " << momentum[p][xHyH + z] << std::endl;
              }
            }
          }
          ++p;
        }
      }
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_momenta\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
} // internal namespace ends here
/******************************************************************************/
/******************************************************************************/
// constructor ****************************************************************/
BasicOperator::BasicOperator () {
  try{
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int dim_row = global_data->get_dim_row();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int number_of_max_mom = global_data->get_number_of_max_mom();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = quarks[0].number_of_dilution_T
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int verbose = global_data->get_verbose();
    // creating gamma matrices
    gamma = new Eigen::SparseMatrix<std::complex<double> >[20];
    for(int i = 0; i < 20; ++i){
      (gamma[i]).resize(4, 4); // TODO: Might not be necessary. CHECK!
      create_gamma(gamma, i);
    }
    if(verbose) {
      std::cout << "gamma 0:\n" << gamma[0] << std::endl;
      std::cout << "gamma 1:\n" << gamma[1] << std::endl;
      std::cout << "gamma 2:\n" << gamma[2] << std::endl;
      std::cout << "gamma 3:\n" << gamma[3] << std::endl;
    }
    // Initializing memory for eigen vectors
    V = new Eigen::MatrixXcd[Lt];
    for(int t = 0; t < Lt; ++t)
      V[t] = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);
    // momentum creation
    momentum = new std::complex<double>*[number_of_max_mom];
    for(int p = 0; p < number_of_max_mom; ++p)
      momentum[p] = new std::complex<double>[Lx * Ly * Lz];
    create_momenta(momentum);
    // memory for the perambulator, random vector and basic operator
    perambulator = new Eigen::MatrixXcd[number_of_rnd_vec];
    rnd_vec = new Eigen::VectorXcd[number_of_rnd_vec];
    op_tsource_tsink = new Eigen::MatrixXcd[number_of_rnd_vec];
    op_tsink_tsource = new Eigen::MatrixXcd[number_of_rnd_vec];
    basicoperator = new Eigen::MatrixXcd**[number_of_rnd_vec];
    contraction = new Eigen::MatrixXcd*[number_of_rnd_vec];
    contraction2 = new Eigen::MatrixXcd*[number_of_rnd_vec];
for(int i = 0; i < number_of_rnd_vec; ++i){
      perambulator[i] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec * Lt,
          number_of_inversions);
      rnd_vec[i] = Eigen::VectorXcd::Zero(Lt * number_of_eigen_vec * 4);
      op_tsink_tsource[i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
      op_tsource_tsink[i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
      basicoperator[i] = new Eigen::MatrixXcd*[Lt];
      contraction[i] = new Eigen::MatrixXcd[16];
      contraction2[i] = new Eigen::MatrixXcd[16];
     for(int t = 0; t < Lt; ++t){
        basicoperator[i][t] = new Eigen::MatrixXcd[4];
        for(int blocknr = 0; blocknr < 4; ++blocknr){
          basicoperator[i][t][blocknr] =  Eigen::MatrixXcd::Zero(
              quarks[0].number_of_dilution_E,
              number_of_eigen_vec);
          // blocks in Basicoperator are on diagonal in the beginning. 
          // index = col + 4 * row
        }
        for(int dirac = 0; dirac < 16; ++dirac){
          contraction[i][dirac] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec, number_of_eigen_vec);
          contraction2[i][dirac] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec, number_of_eigen_vec);
       }
      }
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::BasicOperator\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
// destructor *****************************************************************/
BasicOperator::~BasicOperator () {

  try{
    delete[] perambulator;
    delete[] diluted_operator;
    delete[] gamma;
    delete[] V;

    gamma = NULL;
    V = NULL;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::~BasicOperator\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::read_eigenvectors_from_file (const int config_i) {

  try{
    clock_t t = clock();
    const int Lt = global_data->get_Lt();
    const int dim_row = global_data->get_dim_row();
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();

    std::string filename = global_data->get_path_eigenvectors() + "/"
        + global_data->get_name_eigenvectors();
    //buffer for read in
    std::complex<double>* eigen_vec = new std::complex<double>[dim_row];

    if(verbose) printf("reading eigen vectors from files:\n");
    else printf("\treading eigenvectors:");
    fflush(stdout);

    for(int t = 0; t < Lt; ++t){
      //setting up file
      char name[200];
      sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t); 
      if(verbose) std::cout << "Reading file: " << name << std::endl;
      std::ifstream infile(name, std::ifstream::binary);
    
      for (int nev = 0; nev < number_of_eigen_vec; ++nev) {
        infile.read( (char*) eigen_vec, 2*dim_row*sizeof(double));
        for(int nrow = 0; nrow < dim_row; ++nrow)
        (V[t])(nrow, nev) = eigen_vec[nrow];
      }
      infile.close();
      // small test of trace and sum over the eigen vector matrix!
      if(verbose){
        std::cout << "trace of V^d*V on t = " << t << ":\t"
            << (V[t].adjoint() * V[t]).trace() << std::endl;
        std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"
            << (V[t].adjoint() * V[t]).sum() << std::endl;
      }   
    }
    delete[] eigen_vec;
    t = clock() - t;
    if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::read_eigenvectors_from_file\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::read_perambulators_from_file (const int config_i) {

  try{
    clock_t t = clock();
    char infile[400];
    FILE *fp = NULL;
    const int Lt = global_data->get_Lt();
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = quarks[0].number_of_dilution_T
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int size_perambulator_entry = number_of_inversions * Lt * 4
        * number_of_eigen_vec;

    // memory for reading perambulators
    std::complex<double>* perambulator_read =
        new std::complex<double>[size_perambulator_entry];

    // data path
    std::string filename = global_data->get_path_perambulators() + "/"
        + global_data->get_name_perambulators();

    if(verbose) printf("reading perambulators from files:\n");
    else printf("\treading perambulators:");
    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      //TODO: name is hard-coded at the moment
      sprintf(infile,
          "%s.dil%02d.u.Tso%03d.Dso%01d.Vso%03d.Tsi%03d.Dsi%01d.Vsi%03d.%04d",
          filename.c_str(), rnd_vec_i, quarks[0].number_of_dilution_T,
          quarks[0].number_of_dilution_D, quarks[0].number_of_dilution_E, Lt, 4,
          number_of_eigen_vec, config_i);
      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      if(verbose) printf("\tread file: %s\n", infile);

//      read_lime_spinor((double*) perambulator_read, infile, 0,
//          number_of_inversions * Lt, 2 * Lt * number_of_eigen_vec * 4);
      //TODO: MUST BE CHANGED TO LIME STUFF!!
      fread(perambulator_read, sizeof(std::complex<double>),
          size_perambulator_entry, fp);

      // copy into matrix structure
      int col_i, row_i;
#pragma omp parallel for private(col_i, row_i) schedule(guided)
//      for(row_i = 0; row_i < 4 * Lt * number_of_eigen_vec; ++row_i)
      for(int t1 = 0; t1 < Lt; ++t1)
        for(int ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(int dirac1 = 0; dirac1 < 4; ++dirac1)
//            for(col_i = 0; col_i < number_of_inversions; ++col_i){
            for(int t2 = 0; t2 < quarks[0].number_of_dilution_T; ++t2)
              for(int ev2 = 0; ev2 < quarks[0].number_of_dilution_E; ++ev2)
                for(int dirac2 = 0; dirac2 < quarks[0].number_of_dilution_D; ++dirac2){
              row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
              col_i = quarks[0].number_of_dilution_D * quarks[0].number_of_dilution_E * t2 + quarks[0].number_of_dilution_D * ev2 + dirac2;
              perambulator[rnd_vec_i](4 * number_of_eigen_vec * t1 + number_of_eigen_vec * dirac1 + ev1, 
                quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D * t2 + quarks[0].number_of_dilution_E * dirac2 + ev2) = 
                perambulator_read[row_i * number_of_inversions + col_i];

            }
      fclose(fp);
    }
    delete[] perambulator_read;
    t = clock() - t;
    if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }
  catch(std::exception& e){
    std::cout << e.what()
        << "in: BasicOperator::read_perambulators_from_file\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::read_rnd_vectors_from_file (const int config_i) {

  try{
    clock_t t = clock();
    char infile[400];
    FILE *fp = NULL;
    const int Lt = global_data->get_Lt();
    const int verbose = global_data->get_verbose();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int rnd_vec_length = Lt * number_of_eigen_vec * 4;
    // memory for reading random vectors
    std::complex<double>* rnd_vec_read =
        new std::complex<double>[rnd_vec_length];
    std::string filename = global_data->get_path_perambulators()
        + "/randomvector";

    if(verbose) printf("reading random vectors from files:\n");
    else printf("\treading random vectors:");
    int check_read_in = 0;

    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      // read random vector
      sprintf(infile, "%s.%03d.u.Ti.%04d", filename.c_str(), rnd_vec_i,
          config_i);
      if(verbose) printf("\tread file: %s\n", infile);
      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      check_read_in += fread(rnd_vec_read, sizeof(std::complex<double>),
          rnd_vec_length, fp);
      // copy into matrix structure
      for(int row_i = 0; row_i < rnd_vec_length; ++row_i){
        rnd_vec[rnd_vec_i](row_i) = rnd_vec_read[row_i];
      }
      fclose(fp);
    }
    delete[] rnd_vec_read;
    t = clock() - t;
    if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::read_eigenvectors_from_file\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::build_source_matrix () {

    clock_t t = clock();
    printf("\tbuild source matrix:");
    fflush(stdout);

    const int Lt = global_data->get_Lt();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int dim_row = global_data->get_dim_row();
    const int Vs = global_data->get_Lx() * global_data->get_Ly()
        * global_data->get_Lz();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

    // TODO: checking the order of loops - enhancement might be possible
    for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
      for(int t = 0; t < Lt; ++t){
        // intermediate memory
        // source is a matrix where the source on one timeslice is stored
//        Eigen::MatrixXcd source = Eigen::MatrixXcd::Zero(4 * dim_row,
//            quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D);
//        // V_mat is the EV-matrix with Dirac components: diagonal in Dirac space
//          Eigen::MatrixXcd V_mat = Eigen::MatrixXcd::Zero(4 * dim_row,
//              4 * number_of_eigen_vec);
//        // building matrix of sources on time slice t
//      for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i){
//        for(int mu = 0; mu < 4; ++mu){
//          // filling V_vec and V_mat for source and operator creation
//            Eigen::VectorXcd V_vec = Eigen::VectorXcd::Zero(4 * dim_row);
//// *******************************************************************
//#pragma omp parallel for
//            for(int xs = 0; xs < Vs; ++xs){
//              // three colour components
//              V_vec[xs * 12 + mu * 3 + 0] = (V[t])(xs * 3 + 0, vec_i);
//              V_vec[xs * 12 + mu * 3 + 1] = (V[t])(xs * 3 + 1, vec_i);
//              V_vec[xs * 12 + mu * 3 + 2] = (V[t])(xs * 3 + 2, vec_i);
//            }
//// *******************************************************************
//            V_mat.col(4 * vec_i + mu) = V_vec;
//            // filling the source vector
//            source.col(mu + 4 * (vec_i % quarks[0].number_of_dilution_E)) +=
//                V_vec
//                    * rnd_vec[rnd_i](
//                        mu + vec_i * 4 + 4 * number_of_eigen_vec * t);
//
//
//          }
//        }
        // creating basic operator

      for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i)
        for(int blocknr = 0; blocknr < 4; ++blocknr)
          // blocknr is identical with dirac. basicoperator blockdiagonal in diracspace
          // thus, one can "sort" by dirac index
          ((basicoperator[rnd_i][t][blocknr]))(
                  vec_i % quarks[0].number_of_dilution_E, 
                  vec_i) = 
                  std::conj(rnd_vec[rnd_i](blocknr + vec_i * 4 + 4 * number_of_eigen_vec * t));

//        (basicoperator[rnd_i][t][4]).noalias() = source.adjoint() * V_mat;
      } // loop over time ends here
    }
    t = clock() - t;
    printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }


void BasicOperator::init_operator (const int t_source, const int t_sink){

  clock_t t = clock();
  int i, j;
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i)
    for(int dirac = 0; dirac < 16; ++dirac)
 
  
        // TODO: It might be necessary to change the parrallel part in 'eigen'. CHECK!
        //#pragma omp parallel for shared(op_D) \
        //  private(row_i, col_i) schedule(guided)
      for(int row = 0; row  < 4; ++row){
        i = dirac / 4;
        j = dirac % 4;
        // right multiplication by dirac because [V,dirac] = 0
        contraction[rnd_i][dirac].block(row * number_of_eigen_vec, 0, number_of_eigen_vec, number_of_eigen_vec) =
          perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + number_of_eigen_vec * row,
            quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D * (t_sink / quarks[0].number_of_dilution_T) + quarks[0].number_of_dilution_E * i,
            number_of_eigen_vec,
            quarks[0].number_of_dilution_E) *
          basicoperator[rnd_i][t_sink][j];
        contraction2[rnd_i][dirac].block(row * number_of_eigen_vec, 0, number_of_eigen_vec, number_of_eigen_vec) =
          (perambulator[rnd_i].block(4 * number_of_eigen_vec * t_source + number_of_eigen_vec * i,
            quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D * (t_sink / quarks[0].number_of_dilution_T) + quarks[0].number_of_dilution_E * row,
            number_of_eigen_vec,
            quarks[0].number_of_dilution_E) *
          basicoperator[rnd_i][t_sink][row]).adjoint();
      }

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  }

void BasicOperator::get_operator (){
  clock_t t = clock();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){

 // Eigen::MatrixXcd op_D = Eigen::MatrixXcd::Zero(
 //   4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
  
  int dirac = 5;

  switch(dirac) {
    
    case 5:
      (op_tsource_tsink[rnd_i]).block(0, 0,                       4 * number_of_eigen_vec, number_of_eigen_vec) = contraction[rnd_i][0];
      (op_tsource_tsink[rnd_i]).block(0, 1 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = contraction[rnd_i][5];
      (op_tsource_tsink[rnd_i]).block(0, 2 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = -1 * contraction[rnd_i][10];
      (op_tsource_tsink[rnd_i]).block(0, 3 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = -1 * contraction[rnd_i][15];

      (op_tsink_tsource[rnd_i]).block(0, 0                      , 4 * number_of_eigen_vec, number_of_eigen_vec) = contraction2[rnd_i][0];
      (op_tsink_tsource[rnd_i]).block(0, 1 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = contraction2[rnd_i][5];
      (op_tsink_tsource[rnd_i]).block(0, 2 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = -1 * contraction2[rnd_i][10];
      (op_tsink_tsource[rnd_i]).block(0, 3 * number_of_eigen_vec, 4 * number_of_eigen_vec, number_of_eigen_vec) = -1 * contraction2[rnd_i][15];
      (op_tsink_tsource[rnd_i]).block(0, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec) *= -1;
      (op_tsink_tsource[rnd_i]).block(2 * number_of_eigen_vec, 0, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec) *= -1;

      break;

    default:
      printf("Gamma structure %d not found\n", dirac);
      exit(0);
    } 

  //gamma_5 trick multiplies 2,3,6,7,8,9,12,13 with -1

//    op_tsink_tsource[rnd_i] = (op_tsource_tsink[rnd_i]);
 //   (op_tsink_tsource[rnd_i]).block(0, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec) *= -1;
 //   (op_tsink_tsource[rnd_i]).block(2 * number_of_eigen_vec, 0, 2 * number_of_eigen_vec, 2 * number_of_eigen_vec) *= -1;

  t = clock() - t;
  //printf("\t\tSUCCESS - %.1f seconds\n", ((float) t)/CLOCKS_PER_SEC);
  
  }
  }

  
 

#if 0
#pragma omp parallel for
      for(int t = 0; t < Lt; ++t)
        for(int dirac = 0; dirac < 16; ++dirac)
          if(dirac != 4) basicoperator[rnd_i][t][dirac] = create_operator(t,
              rnd_i, dirac, 0);
    } // loop over random vectors ends here
   /******************************************************************************/
/******************************************************************************/
/******************************************************************************/
Eigen::MatrixXcd BasicOperator::mul_l_gamma (const Eigen::MatrixXcd& in,
    const int rows, const int cols, const int mu, const int flag) const {

    int row, col;
    Eigen::SparseMatrix<std::complex<double> > dirac;
    Eigen::MatrixXcd op_D = Eigen::MatrixXcd::Zero(rows, cols);

    switch(flag) { // **********************************************************
    case 0:
      dirac = gamma[mu];
      break;
    case 1:
      dirac = gamma[0] * gamma[mu].adjoint() * gamma[0];
      break;
    case 2:
      dirac = gamma[0] * gamma[mu] * gamma[0];
      break;
    case 3:
      dirac = gamma[5] * gamma[mu];
      break;
    case 4:
      dirac = gamma[5] * gamma[mu] * gamma[5];
      break;
    default:
      printf("Operator construction flag not found\n");
      exit(0);
    } // switch ends here ******************************************************

// TODO: It might be necessary to change the parrallel part in 'eigen'. CHECK!
#pragma omp parallel for shared(op_D, dirac) \
    private(row, col) schedule(guided)
    for(row = 0; row < rows / 4; ++row)
      for(col = 0; col < cols / 4; ++col)
        // right multiplication by dirac because [V,dirac] = 0
        op_D.block(4 * row, 4 * col, 4, 4) = dirac
            * in.block(4 * row, 4 * col, 4, 4);

    return op_D;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
Eigen::MatrixXcd BasicOperator::mul_r_gamma (const Eigen::MatrixXcd& in,
    const int rows, const int cols, const int mu, const int flag) const {

    int row, col;
    Eigen::SparseMatrix<std::complex<double> > dirac1;
    Eigen::MatrixXcd dirac = Eigen::MatrixXcd::Zero(cols, cols);
    Eigen::Matrix4cd temp;
    Eigen::MatrixXcd op_D = Eigen::MatrixXcd::Zero(rows, cols);

    switch(flag) { // **********************************************************
    case 0:
      dirac1 = gamma[mu];
      break;
    case 1:
      dirac1 = gamma[0] * gamma[mu].adjoint() * gamma[0];
      break;
    case 2:
      dirac1 = gamma[0] * gamma[mu] * gamma[0];
      break;
    case 3:
      dirac1 = gamma[5] * gamma[mu];
      break;
    case 4:
      dirac1 = gamma[5] * gamma[mu] * gamma[5];
      break;
    default:
      printf("Operator construction flag not found\n");
      exit(0);
    } // switch ends here ******************************************************

// TODO: It might be necessary to change the parrallel part in 'eigen'. CHECK!

    int r = rows / 4;
    int c = cols / 4;

    temp = dirac1;

    for(row = 0; row < cols; ++row)
      for(col = 0; col < cols; ++col)
        dirac(row, cols) = temp(row / r, col / c);

#pragma omp parallel for shared(op_D, dirac) \
    private(row, col) schedule(guided)
    for(row = 0; row < 4; ++row){
      for(col = 0; col < 4; ++col){
        // right multiplication by dirac because [V,dirac] = 0
        op_D = in * dirac;
        }
    }

    return op_D;
}

#endif
