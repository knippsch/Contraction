//everything to read and write from/to files

#include "ReadWrite.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
// constructor ****************************************************************/
/******************************************************************************/
/******************************************************************************/
ReadWrite::ReadWrite () : perambulator(), 
                          rnd_vec() {

  const int Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const int number_of_inversions = (Lt / quarks[0].number_of_dilution_T)
      * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;

  // memory for perambulator and random vector
  perambulator.resize(number_of_rnd_vec);
  for(int i = 0; i < number_of_rnd_vec; ++i){
    perambulator[i] = Eigen::MatrixXcd::Zero(4 * number_of_eigen_vec * Lt,
        number_of_inversions);
  }
  // TODO: the resize is unnasassarry if it is done in initialiser list 
  // with the ctor
  rnd_vec.resize(number_of_rnd_vec, 
                 LapH::RandomVector(Lt*number_of_eigen_vec*4));
    
  std::cout << "\tallocated memory for ReadWrite" << std::endl;

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
      sprintf(temp, "cnfg%d/rnd_vec_%01d/", config_i, rnd_vec_i);
      std::string filename = global_data->get_path_perambulators() + "/"
          + temp;

      // data path for juqueen contractions
//      sprintf(temp, "cnfg%d/", config_i);
//      std::string filename = global_data->get_path_perambulators() + "/"
//          + temp;

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
      int blabla = fread(&(perambulator_read[0]), sizeof(std::complex<double>),
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

  clock_t t = clock();
  char infile[400];
  const int Lt = global_data->get_Lt();
  const int verbose = global_data->get_verbose();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const int rnd_vec_length = Lt * number_of_eigen_vec * 4;

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
//				+ "/" + temp;

    // read random vector
    sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev%04d.%04d", 
        filename.c_str(), rnd_vec_i, number_of_eigen_vec, config_i);
//      sprintf(infile, "%srandomvector.rndvecnb%02d.u.nbev0120.%04d", 
//          filename.c_str(), rnd_vec_i, config_i);

//      sprintf(infile, "%s.%03d.u.Ti.%04d", filename.c_str(), rnd_vec_i,
//          config_i);

    // TODO:: explicit type conversion - Bad style
    rnd_vec[rnd_vec_i].read_random_vector(infile);
  }
  t = clock() - t;
  if(!verbose) std::cout << "\t\tSUCCESS - " << std::fixed 
                         << std::setprecision(1)
                         << ((float) t)/CLOCKS_PER_SEC << " seconds" 
                         << std::endl; 
}
