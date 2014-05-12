//everything to read and write from/to files

#include "ReadWrite.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);

static void create_momenta (std::complex<double>** momentum) {

  try{
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int number_of_max_mom = global_data->get_number_of_max_mom();
    const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();
    // helper variables for momenta
    const double px = 2. * M_PI / (double) Lx;
    const double py = 2. * M_PI / (double) Ly;
    const double pz = 2. * M_PI / (double) Lz;
    int p = 0;
    int max_mom_squared = number_of_max_mom * number_of_max_mom;

    // running over all momentum components
    for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
      for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
        for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
          if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
            continue;
          }

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
    std::cout << e.what() << "in: ReadWrite::create_momenta\n";
    exit(0);
  }
}

int check_momenta() {
  try {

    const int number_of_max_mom = global_data->get_number_of_max_mom();
    const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();

    int p = 0;
    int max_mom_squared = number_of_max_mom * number_of_max_mom;

    // running over all momentum components
    for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
      for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
        for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
          if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
            continue;
          }
          p++;
        }
      }
    }

    return p;
  }
  catch(std::exception& e) {
    std::cout << e.what() << "in: ReadWrite::check_momenta\n";
    exit(0);
  }
}
   

/******************************************************************************/
/******************************************************************************/
} // internal namespace ends here

/******************************************************************************/
/******************************************************************************/
// constructor ****************************************************************/
/******************************************************************************/
/******************************************************************************/

ReadWrite::ReadWrite () {
  try{
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int dim_row = global_data->get_dim_row();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_max_mom = global_data->get_number_of_max_mom();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = (Lt / quarks[0].number_of_dilution_T)
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;

    // Initializing memory for eigen vectors
    V = new Eigen::MatrixXcd[Lt];
    for(int t = 0; t < Lt; ++t)
      V[t] = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);
    V_temp = new Eigen::MatrixXcd[Lt];
    for(int t = 0; t < Lt; ++t)
      V_temp[t] = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);

    // momentum creation
    number_of_momenta = check_momenta();
    std::cout << "Calulatiing correlators for " << number_of_momenta << 
    " momenta" << std::endl;
    momentum = new std::complex<double>*[number_of_momenta];
    for(int p = 0; p < number_of_momenta; ++p)
      momentum[p] = new std::complex<double>[Lx * Ly * Lz];
    create_momenta(momentum);

    // memory for the perambulator, random vector and basic operator
    basicoperator = new Eigen::MatrixXcd***[number_of_momenta];
    for(int p = 0; p < number_of_momenta; ++p) {
      basicoperator[p] = new Eigen::MatrixXcd**[number_of_rnd_vec];
      for(int i = 0; i < number_of_rnd_vec; ++i){
        basicoperator[p][i] = new Eigen::MatrixXcd*[Lt];
        for(int t = 0; t < Lt; ++t){
          basicoperator[p][i][t] = new Eigen::MatrixXcd[4];
          for(int blocknr = 0; blocknr < 4; ++blocknr){
            // blocks in Basicoperator are on diagonal in the beginning. 
            // non-zero blocks have row = col = blocknr
            basicoperator[p][i][t][blocknr] =  Eigen::MatrixXcd::Zero(
                quarks[0].number_of_dilution_E,
                number_of_eigen_vec);
          }
        }
      }   
    }

    perambulator = new Eigen::MatrixXcd[number_of_rnd_vec];
    rnd_vec = new Eigen::VectorXcd[number_of_rnd_vec];
    for(int i = 0; i < number_of_rnd_vec; ++i){
      perambulator[i] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec * Lt,
          number_of_inversions);
      rnd_vec[i] = Eigen::VectorXcd::Zero(Lt * number_of_eigen_vec * 4);
    }
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::ReadWrite\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
// destructor *****************************************************************/
/******************************************************************************/
/******************************************************************************/

ReadWrite::~ReadWrite() {

  try{
    delete[] perambulator;
    delete[] basicoperator;
    delete[] rnd_vec;
    delete[] V;
    delete[] momentum;

    V = NULL;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::~ReadWrite\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::build_source_matrix (const int p) {

  clock_t t2 = clock();
  printf("\tbuild source matrix:\n");
  fflush(stdout);

  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int dim_row = global_data->get_dim_row();
  //const int Vs = global_data->get_Lx() * global_data->get_Ly()
  //    * global_data->get_Lz();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const int number_of_max_mom = global_data->get_number_of_max_mom();

  // creating basic operator

  // for p = 0, s is the unit matrix. Thus, the V.adjoint() * V multiplication
  // can be omitted
  if(p == 0) {
    for(int t = 0; t < Lt; ++t){
      for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i) {
        for(int blocknr = 0; blocknr < 4; ++blocknr) {
          for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
            // blocknr is identical with dirac. basicoperator blockdiagonal 
            // in diracspace -> treat every dirac index individually
            ((basicoperator[0][rnd_i][t][blocknr]))(
                vec_i % quarks[0].number_of_dilution_E, vec_i) =
                std::conj(rnd_vec[rnd_i](blocknr + vec_i * 4 + 4 * 
                number_of_eigen_vec * t));
          }
        }
      }
    } // loop over time ends here
  }

  // case p != 0
  else  {
    // TODO: checking the order of loops - enhancement might be possible
    for(int t = 0; t < Lt; ++t){
    
      // multiply EV with momenta
      for(int x = 0; x < dim_row; ++x) {
        V_temp[t].row(x) = momentum[p][x/3] * V[t].row(x);
      }
      
      //TODO: check if saving V[t].adjoint in own matrix is faster
      Eigen::MatrixXcd s = (V[t]).adjoint() * V_temp[t];

      for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
        for(int blocknr = 0; blocknr < 4; ++blocknr) {
          // blocknr is identical with dirac. basicoperator blockdiagonal in 
          // diracspace -> treat every dirac index individually
          for(int vec_i = 0; vec_i < number_of_eigen_vec; ++vec_i) {
            basicoperator[p][rnd_i][t][blocknr].row(vec_i % 
                quarks[0].number_of_dilution_E) +=
                std::conj(rnd_vec[rnd_i](blocknr + vec_i * 4 + 
                4 * number_of_eigen_vec * t)) * s.row(vec_i);
          }
        }
      }
    } // loop over time ends here
  }

  t2 = clock() - t2;
  printf("\t\tSUCCESS - %.1f seconds\n", ((float) t2)/CLOCKS_PER_SEC);
  fflush(stdout);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_eigenvectors_from_file (const int config_i) {

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
    std::cout << e.what() << "in: ReadWrite::read_eigenvectors_from_file\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_perambulators_from_file (const int config_i) {

  try{
    clock_t t = clock();
    char infile[400];
    FILE *fp = NULL;
    const int Lt = global_data->get_Lt();
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = (Lt / quarks[0].number_of_dilution_T)
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int size_perambulator_entry = number_of_inversions * Lt * 4
        * number_of_eigen_vec;

    // memory for reading perambulators
    std::complex<double>* perambulator_read =
        new std::complex<double>[size_perambulator_entry];
    char temp[9];

    if(verbose){
      printf("reading perambulators from files:\n");
    }
    else{
      printf("\treading perambulators:");
    }

    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      // data path
      sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
      std::string filename = global_data->get_path_perambulators() + "/"
          + temp;

      //TODO: name is hard-coded at the moment
      sprintf(infile,
          "%sperambulator.rndvecnb%02d.u.TsoB0024.VsoI0006.DsoF4.TsiF0048."
          "SsiF13824.DsiF4.CsiF3.smeared1.%05d", filename.c_str(), rnd_vec_i, 
          config_i);

//         "%s.dil%02d.u.Tso%03d.Dso%01d.Vso%03d.Tsi%03d.Dsi%01d.Vsi%03d.%04d",
//          filename.c_str(), rnd_vec_i, quarks[0].number_of_dilution_T,
//          quarks[0].number_of_dilution_D, quarks[0].number_of_dilution_E, Lt, 4,
//          number_of_eigen_vec, config_i);

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
      for(int t1 = 0; t1 < Lt; ++t1)
        for(int ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(int dirac1 = 0; dirac1 < 4; ++dirac1)
            for(int t2 = 0; t2 < (Lt / quarks[0].number_of_dilution_T); ++t2)
              for(int ev2 = 0; ev2 < quarks[0].number_of_dilution_E; ++ev2)
                for(int dirac2 = 0; dirac2 < quarks[0].number_of_dilution_D; 
                    ++dirac2){
                  row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
                  col_i = quarks[0].number_of_dilution_D * 
                      quarks[0].number_of_dilution_E * t2 + 
                      quarks[0].number_of_dilution_D * ev2 + dirac2;
                  perambulator[rnd_vec_i](4 * number_of_eigen_vec * t1 + 
                      number_of_eigen_vec * dirac1 + ev1, 
                      quarks[0].number_of_dilution_E * 
                        quarks[0].number_of_dilution_D * t2 + 
                        quarks[0].number_of_dilution_E * dirac2 + ev2) = 
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
        << "in: ReadWrite::read_perambulators_from_file\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_rnd_vectors_from_file (const int config_i) {

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
    char temp[9];

    if(verbose){
      printf("reading random vectors from files:\n");
    }
    else{
      printf("\treading random vectors:");
    }

    int check_read_in = 0;

    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      // data path
      sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
      std::string filename = global_data->get_path_perambulators()
				+ "/" + temp;

      // read random vector
      sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev0120.%04d", 
          filename.c_str(), rnd_vec_i, config_i);

//      sprintf(infile, "%s.%03d.u.Ti.%04d", filename.c_str(), rnd_vec_i,
//          config_i);

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
    std::cout << e.what() << "in: ReadWrite::read_eigenvectors_from_file\n";
    exit(0);
  }
}
