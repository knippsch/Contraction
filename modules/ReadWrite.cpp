//everything to read and write from/to files

#include "ReadWrite.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);

static void create_momenta (array_cd_d2 momentum) {

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
          //TODO: for Lx == Ly == Lz ipxH and ipxHipyH may be integers and px, 
          //py get multiplied in the exponential
          // running over all lattice points
          for(int x = 0; x < Lx; ++x){
            const int xH = x * Ly * Lz; // helper variable
            const double ipxH = ipx * px * x; // helper variable
            for(int y = 0; y < Ly; ++y){
              const int xHyH = xH + y * Lz; // helper variable
              const double ipxHipyH = ipxH + ipy * py * y; // helper variable
              for(int z = 0; z < Lz; ++z){
                momentum[p][xHyH + z] = exp(-I * (ipxHipyH + ipz * pz * z));
              }
            }
          }
          ++p;
        }
      }
    }
  return;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::create_momenta\n";
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

ReadWrite::ReadWrite () : perambulator(), rnd_vec(), basicoperator() {
  try{
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int Vs = Lx * Ly * Lz;
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = (Lt / quarks[0].number_of_dilution_T)
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int number_of_momenta = global_data->get_number_of_momenta();
    const std::vector<int> mom_squared = global_data->get_momentum_squared();
    const int V_for_lime = global_data->get_V_for_lime();

    // memory for momentum matrices
    momentum.resize(boost::extents[number_of_momenta][Vs]);
    create_momenta(momentum);

    // memory for basic operator
    basicoperator.resize(boost::extents[number_of_momenta][Lt][4]);
    for(int p = 0; p < number_of_momenta; ++p) {
      for(int t = 0; t < Lt; ++t){
        // changed to case of no displacement. Else dir < 4
        for(int dir = 0; dir < 4; dir++) {
          // blocks in Basicoperator are on diagonal in the beginning. 
          // non-zero blocks have row = col = blocknr
          basicoperator[p][t][dir] =  Eigen::MatrixXcd::Zero(
              number_of_eigen_vec, number_of_eigen_vec);
        }
      }
    }   

    // memory for perambulator and random vector
    perambulator.resize(number_of_rnd_vec);
    rnd_vec.resize(number_of_rnd_vec);
    for(int i = 0; i < number_of_rnd_vec; ++i){
      perambulator[i] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec * Lt,
          number_of_inversions);
      rnd_vec[i] = Eigen::VectorXcd::Zero(Lt * number_of_eigen_vec * 4);
    }

#if 1
    //memory for gauge fields
    gaugefield = new double[V_for_lime];
    //Allocate Eigen Array to hold timeslice
    eigen_timeslice = new Eigen::Matrix3cd *[Vs];
    iup = new int*[Vs];
    idown = new int*[Vs];
    for (int i = 0; i < Vs; ++i ) {
      eigen_timeslice[i] = new Eigen::Matrix3cd[3];
      iup[i] = new int[3];
      idown[i] = new int[3];
    }

    hopping3d(iup, idown);
#endif

  std::cout << "\tallocated memory for ReadWrite" << std::endl;

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

//ReadWrite::~ReadWrite() {
//  try{
//    const int Lt = global_data->get_Lt();
//    const int Lx = global_data->get_Lx();
//    const int Ly = global_data->get_Ly();
//    const int Lz = global_data->get_Lz();
//    const int Vs = Lx * Ly * Lz;
//    const int number_of_momenta = global_data->get_number_of_momenta();
//    // delete all memory
//    delete [] perambulator;
//    delete [] rnd_vec;
//    for(int p = 0; p < number_of_momenta; ++p) {
//      for(int t = 0; t < Lt; ++t){
//        delete [] basicoperator[p][t];
//      }
//      delete [] basicoperator[p];
//    }   
//    delete [] basicoperator;
//    for(int p = 0; p < number_of_momenta; ++p)
//      delete [] momentum[p];
//    delete [] momentum;
//    for (int i = 0; i < Vs; ++i ) {
//      delete [] eigen_timeslice[i];
//      delete [] iup[i];
//      delete [] idown[i];
//    }
//    delete [] eigen_timeslice;
//    delete [] iup;
//    delete [] idown;
//    delete [] gaugefield;
//  }
//  catch(std::exception& e){
//    std::cout << e.what() << "in: ReadWrite::~ReadWrite\n";
//    exit(0);
//  }
//}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


void ReadWrite::build_source_matrix (const int config_i) {

  clock_t t2 = clock();
  std::cout << "\tbuild source matrix:";

  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int dim_row = global_data->get_dim_row();
  const int displ_min = global_data->get_displ_min();
  const int displ_max = global_data->get_displ_max();
  const int number_of_momenta = global_data->get_number_of_momenta();

  // creating basic operator
  Eigen::MatrixXcd V_t = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);
  Eigen::MatrixXcd W_t = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);

  for(int t = 0; t < Lt; ++t){

    (V_t).setZero();
    read_eigenvectors_from_file(V_t, config_i, t);
   
//    std::cout << "V_t with t = " << t << std::endl;
//    std::cout << V_t.block(0,0,6,6) << std::endl;
//    std::cout << basicoperator[number_of_momenta - p - 1][t][0].block(0,0,6,6) << std::endl;
//        << std::endl << "\n" << s.block(0,0,6,6) << std::endl;
//    std::cout << std::endl;

    for(int dir = displ_min; dir < displ_max + 1; dir++) {
      for(int p = number_of_momenta/2; p < number_of_momenta; p++){
        // TODO: implement switch case for displacement
        // case no displacement
        if(dir == 0) {
          // TODO: checking the order of loops - enhancement might be possible
          // e.g. by reordering basicoperator with t faster than rnd_i
      
          // for p = 0, s is the unit matrix. Thus, the V.adjoint() * V multiplication
          // can be omitted
          // TODO: initialize somewhere in the constructor
          if(p == (number_of_momenta/2)) { // zero momentum
            for(int vec_i = 0; vec_i < number_of_eigen_vec; vec_i++) {
              for(int vec_j = 0; vec_j < number_of_eigen_vec; vec_j++) {
                (basicoperator[number_of_momenta/2][t][0])(vec_i, vec_i) = 1;
              }
            }
          } else { // not zero momentum
            // momentum vector contains exp(-i p x)
            // Divisor 3 for colour index. All three colours on same lattice site get
            // the same momentum
            for(int x = 0; x < dim_row; ++x) {
              mom(x) = momentum[p][x/3];
            }
            basicoperator[p][t][0] = V_t.adjoint() * mom.asDiagonal() * V_t;
            basicoperator[number_of_momenta - p - 1][t][0] = 
                (basicoperator[p][t][0]).adjoint();
          } // end if momentum
        // case displacement
        } else {
      
          (W_t).setZero();

          //Time Slice of Configuration
          //Factor 3 for second color index, 2 for complex numbers
          double* timeslice = gaugefield + (t * dim_row * 3 * 2);
            
          //Write Timeslice in Eigen Array
          map_timeslice_to_eigen(eigen_timeslice, timeslice);
        
          //displacement in one direction i acting to the right
          right_displacement_one_dir(eigen_timeslice, iup, idown, dir - 1, 
              V_t, W_t);
          // dir = 3 for z-displacement (1 x 2 y)
          // W holds DV
  
          if(p == number_of_momenta/2) {
            
            basicoperator[number_of_momenta/2][t][dir] = V_t.adjoint() * W_t
                - W_t.adjoint() * V_t;

          } else {   

            // momentum vector contains exp(-i p x)
            // Divisor 3 for colour index. All three colours on same lattice site get
            // the same momentum
            for(int x = 0; x < dim_row; ++x) {
              mom(x) = momentum[p][x/3];
            }
              
            //TODO: check if saving V[t].adjoint in own matrix is faster
            //TODO: is that efficient?
            // build basicoperator = V^dagger exp(-ipx) V 
            // opposite momentum is just basicoperator daggered
            basicoperator[p][t][dir] = V_t.adjoint() * mom.asDiagonal() *  W_t - 
                W_t.adjoint() * mom.asDiagonal() * V_t;
            basicoperator[number_of_momenta - p - 1][t][dir] = 
                (-1) *  (basicoperator[p][t][dir]).adjoint();
          
          }
        } // end if displacement
  
      } // end for momentum
    } // end for displacement

  } // loop over time ends here

  t2 = clock() - t2;
  std::cout << "\t\tSUCCESS - " << std::fixed << std::setprecision(1)
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
  return;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_eigenvectors_from_file (Eigen::MatrixXcd& V, const int config_i, const int t) {

  try{
    //clock_t time = clock();
 //   const int Lt = global_data->get_Lt();
    const int dim_row = global_data->get_dim_row();
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();

    std::string filename = global_data->get_path_eigenvectors() + "/"
        + global_data->get_name_eigenvectors();
    //buffer for read in
    vec eigen_vec(dim_row);

    //if(verbose) printf("reading eigen vectors from files:\n");
    //else printf("\treading eigenvectors:");
    //fflush(stdout);

    //setting up file
    char name[200];
    sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t); 
    if(verbose) std::cout << "Reading file: " << name << std::endl;
    std::ifstream infile(name, std::ifstream::binary);
  
    for (int nev = 0; nev < number_of_eigen_vec; ++nev) {
      infile.read( (char*) &(eigen_vec[0]), 2*dim_row*sizeof(double));
      for(int nrow = 0; nrow < dim_row; ++nrow){
        V(nrow, nev) = eigen_vec[nrow];
      }
    }
    infile.close();
    // small test of trace and sum over the eigen vector matrix!
    if(verbose){
      std::cout << "trace of V^d*V on t = " << t << ":\t"
          << (V.adjoint() * V).trace() << std::endl;
      std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"
          << (V.adjoint() * V).sum() << std::endl;
    }   

    //time = clock() - time;
    //if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
  return;
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
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int Vs = Lx * Ly * Lz;
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_inversions = (Lt / quarks[0].number_of_dilution_T)
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
    const int size_perambulator_entry = number_of_inversions * Lt * 4
        * number_of_eigen_vec;

    // memory for reading perambulators
    vec perambulator_read(size_perambulator_entry);
    char temp[100];

    if(verbose){
      std::cout << "\treading perambulators from files:" << std::endl;
    } else {
      std::cout << "\treading perambulators:";
    }

    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      //data path for christians perambulators
//      std::string filename = global_data->get_path_perambulators() + "/";

      // data path for qbig contractions
//      sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
//      std::string filename = global_data->get_path_perambulators() + "/"
//          + temp;

      // data path for juqueen contractions
      sprintf(temp, "cnfg%d/", config_i);
      std::string filename = global_data->get_path_perambulators() + "/"
          + temp;

      //TODO: sink dilution is  hard-coded at the moment
      sprintf(infile,
          "%sperambulator.rndvecnb%02d.u.TsoB%04d.VsoI%04d.DsoF%1d.TsiF%04d."
          "SsiF%d.DsiF4.CsiF3.smeared1.%05d", 
          filename.c_str(), rnd_vec_i, 
          Lt / quarks[0].number_of_dilution_T, quarks[0].number_of_dilution_E,
          quarks[0].number_of_dilution_D,
          Lt, Vs,
          config_i);

      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      if(verbose) std::cout << "\tread file" << infile << std::endl;

//      read_lime_spinor((double*) perambulator_read, infile, 0,
//          number_of_inversions * Lt, 2 * Lt * number_of_eigen_vec * 4);
      //TODO: MUST BE CHANGED TO LIME STUFF!!
      fread(&(perambulator_read[0]), sizeof(std::complex<double>),
          size_perambulator_entry, fp);

      // copy into matrix structure
      int col_i, row_i;
      int t1, t2, ev1, ev2, dirac1, dirac2;
#pragma omp parallel for private(col_i, row_i, t1, t2, ev1, ev2, dirac1, dirac2) schedule(guided)
      for(t1 = 0; t1 < Lt; ++t1)
        for(ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(dirac1 = 0; dirac1 < 4; ++dirac1)
            for(t2 = 0; t2 < (Lt / quarks[0].number_of_dilution_T); ++t2)
              for(ev2 = 0; ev2 < quarks[0].number_of_dilution_E; ++ev2)
                for(dirac2 = 0; dirac2 < quarks[0].number_of_dilution_D; ++dirac2){
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
    t = clock() - t;
    if(!verbose) std::cout << "\t\tSUCCESS - " << std::fixed << std::setprecision(1)
        << ((float) t)/CLOCKS_PER_SEC << " seconds" << std::endl;
  return;
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
    vec rnd_vec_read(rnd_vec_length);
    char temp[100];

    if(verbose){
      std::cout << "\treading randomvectors from files:" << std::endl;
    } else {
      std::cout << "\treading randomvectors:";
    }

    int check_read_in = 0;
    for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
      // data path Christians perambulators
//      std::string filename = global_data->get_path_perambulators() + "/";

      // data path for qbig contractions
      sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
      std::string filename = global_data->get_path_perambulators()
        + "/" + temp;

      // data path for juqueen contractions
//      sprintf(temp, "cnfg%d/", config_i);
//      std::string filename = global_data->get_path_perambulators()
//        + "/" + temp;

      // data path for juqueen contractions
      sprintf(temp, "cnfg%d/", config_i);
      std::string filename = global_data->get_path_perambulators()
				+ "/" + temp;

      // read random vector
      sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev%04d.%04d", 
          filename.c_str(), rnd_vec_i, number_of_eigen_vec, config_i);

//      sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev0120.%04d", 
//          filename.c_str(), rnd_vec_i, config_i);

//      sprintf(infile, "%s.%03d.u.Ti.%04d", filename.c_str(), rnd_vec_i,
//          config_i);

      if(verbose) printf("\tread file: %s\n", infile);
      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      check_read_in += fread(&(rnd_vec_read[0]), sizeof(std::complex<double>),
          rnd_vec_length, fp);

      // copy into matrix structure
      for(int row_i = 0; row_i < rnd_vec_length; ++row_i){
        rnd_vec[rnd_vec_i](row_i) = rnd_vec_read[row_i];
      }
      fclose(fp);
    }
    t = clock() - t;
    if(!verbose) std::cout << "\t\tSUCCESS - " << std::fixed << std::setprecision(1)
      << ((float) t)/CLOCKS_PER_SEC << " seconds" << std::endl; 
    return;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::read_eigenvectors_from_file\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// Christophers function to read in the gauge fields. Don't ask me for comments on 
// what this does
void ReadWrite::read_lime_gauge_field_doubleprec_timeslices(const int config_i) {

  try{
    clock_t time = clock();
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int verbose = global_data->get_verbose();
    
    const int slice_i = 0;
    const int slice_f = Lt+1;
    char filename[200];

    FILE * ifs;
    int t, x, y, z, status;
    n_uint64_t bytes;
    char * header_type;
    LimeReader * limereader;
    double tmp[72], tmp2[72];
    int words_bigendian;

    if(verbose){
      std::cout << "\treading gauge fields from files:" << std::endl;
    } else {
      std::cout << "\treading gauge fields:";
    }

    sprintf(filename, "%s/conf.%04d",
      global_data->get_config_path().c_str(), config_i);

    words_bigendian = big_endian();
    ifs = fopen(filename, "r");
    if(ifs == (FILE *)NULL) {
      fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
      exit(500);
    }
    limereader = limeCreateReader( ifs );
    if( limereader == (LimeReader *)NULL ) {
      fprintf(stderr, "Unable to open LimeReader\n");
      exit(500);
    }
    if(verbose) std::cout << "\t\tread file: " << filename << std::endl;
    while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
      if(status != LIME_SUCCESS ) {
        fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n",
            status);
        status = LIME_EOF;
        break;
      }
      header_type = limeReaderType(limereader);
      if(strcmp("ildg-binary-data",header_type) == 0) break;
    }
    if(status == LIME_EOF) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
      limeDestroyReader(limereader);
      fclose(ifs);
      exit(-2);
    }
    bytes = limeReaderBytes(limereader);
    if(bytes != (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double)) {
      if(bytes != (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double)/2) {
        fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n",
            (n_uint64_t)bytes, filename,
            (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double));
        fprintf(stderr, "Aborting...!\n");
        fflush( stdout );
        exit(501);
      }
      else {
        fclose(ifs);
        fprintf(stderr, "single precision read!\n");

        fprintf(stderr, "Not implemented!\n");
        exit(EXIT_FAILURE);
        read_lime_gauge_field_singleprec(gaugefield, filename, Lt, Lx, Ly, Lz);
        return;
      }
    }

    bytes = (n_uint64_t)72*sizeof(double);

    for(t = 0; t < Lt; t++) {
      for(z = 0; z < Lz; z++) {
        for(y = 0; y < Ly; y++) {
          for(x = 0; x < Lx; x++) {

            // check for endianess and reading in data
            // the pointer limereader is internally increased by bytes
            // in the limeReaderReadData function
            if(!words_bigendian) {
              status = limeReaderReadData(tmp, &bytes, limereader);
              byte_swap_assign(tmp2, tmp, 72);
            }
            else
              status = limeReaderReadData(tmp2, &bytes, limereader);
            // check if reading was successfull
            if(status < 0 && status != LIME_EOR) {
              fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
                  status, filename);
              exit(500);
            }

            // we just want to read in data in the specific range of timeslices
            // must be here because the limereader pointer must be increased 
            // correctly
            // could be done with much more performance but might be tricky to 
            // do correctly
            if(t<slice_i || t>slice_f)
              continue;

            // copy of link variables from tmp2 into config
            // ILDG has mu-order: x,y,z,t so it is changed here to: t,x,y,z !
            const size_t p = (size_t) ( ((t-slice_i)*Lx*Lz*Lz + 
                x*Ly*Lz + y*Lz + z) * 72); // position in config
            size_t k = 0;
            for(size_t mu = 1; mu <= 4; mu++) { // mu=4 is for the shift of U_t
              size_t index;
              if (mu != 4)
                index = p + mu*18; // for U_x, U_y and U_z
              else
                index = p; // U_t is copied into the beginning of
              // the (config+p) array

              for(size_t i = 0; i < 3; i++) {
                for(size_t j = 0; j < 3; j++) {
                  gaugefield[index+6*i+2*j] = tmp2[2*k];
                  gaugefield[index+6*i+2*j+1] = tmp2[2*k+1];
                  k++;
                }
              }

            } // loop over mu ends here

          } // loop over position space ends here
        }
      }
    }
    if(status < 0 && status != LIME_EOR) {
      fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
          status, filename);
      exit(500);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    time = clock() - time;
    if(!verbose) std::cout << "\t\tSUCCESS - " << std::fixed << std::setprecision(1)
      << ((float) time)/CLOCKS_PER_SEC << " seconds" << std::endl; 
    return;

  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::read_lime_gauge_field_doubleprec_timeslices\n";
    exit(0);
  }

}

