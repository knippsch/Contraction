#include "CorrelatorIo2pt.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//TODO: put that into the IoHelpers. Does not work naively as the template
// type must be known at compile time

// convert multiarray 2pt correlator to vector to match write_2pt_lime
template <typename listvector> 
void convert_C2_mes_to_vec(const listvector& op_mes, array_cd_d2& C2_mes, 
                           std::vector<Tag>& tags, std::vector<vec>& corr){

  const size_t Lt = global_data->get_Lt();

  corr.resize(op_mes.size());
  for(auto& c : corr)
    c.resize(Lt);
  tags.resize(op_mes.size());

  for(const auto& op : op_mes){
    //TODO: solve the copying of C2_mes (boost) into corr (std::vec) without
    // the copy constructor in assign()
    corr[op.id].assign(C2_mes[op.id].origin(), C2_mes[op.id].origin() + 
                       C2_mes[op.id].size());
    set_tag(tags[op.id], op.index.front());

  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void export_corr_2pt(const char* filename, array_cd_d2& C2_mes){
  
  const vec_pdg_C2 op_C2 = global_data->get_op_C2();

  GlobalDat dat;
  std::vector<Tag> tags;
  std::vector<vec> corr;

  convert_C2_mes_to_vec <vec_pdg_C2> (op_C2, C2_mes, tags, corr);
  if (file_exist(filename)){
    char filename_new [150];
    sprintf(filename_new,"_changed");
    std::cout << "file already exists! Renaming to "
    << filename_new << std::endl;
    write_2pt_lime(filename_new, dat, tags, corr);
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void export_corr_4pt(const char* filename, array_cd_d2& C4_mes){
  
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();

  GlobalDat dat;
  std::vector<Tag> tags;
  std::vector<vec> corr;

  convert_C2_mes_to_vec <vec_pdg_C4> (op_C4, C4_mes, tags, corr);
  if (file_exist(filename)){
    char filename_new [150];
    sprintf(filename_new,"_changed");
    std::cout << "file already exists! Renaming to "
    << filename_new << std::endl;
  }
  write_2pt_lime(filename, dat, tags, corr);

}
///////////////////////////////////////////////////////////////////////////////
// Write 2pt correlator ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void write_2pt_lime(const char* filename, GlobalDat& dat, 
                    std::vector<Tag>& tags, std::vector<vec>& corr){

  bool be = big_endian();
  if (be){
    swap_correlators(corr);
    swap_tag_vector(tags);
    dat = swap_glob_dat(dat);
  }
    // calculate checksum of all correlators

    //concatenate all correlation functions in one vector
    std::vector<cmplx> collect;
    for(auto& c : corr)
      for (auto& el : c) collect.push_back(el);
    size_t glob_bytes = collect.size()*2*sizeof(double);
    size_t global_chksum = checksum <std::vector<cmplx> > (collect,
        glob_bytes);
    std::cout << "Global Checksum is: " << global_chksum << std::endl;
    if(be) global_chksum = swap_endian<size_t>(global_chksum);
    write_1st_msg(filename, dat, global_chksum);

  // setup what is needed for output to lime
  FILE* fp;
  LimeWriter* w;
  fp = fopen( filename, "a" );
  w = limeCreateWriter( fp );

  // writing the correlators to the end
  append_msgs(filename, corr, tags, w, be);
  limeDestroyWriter( w );
  fclose(fp);
}
///////////////////////////////////////////////////////////////////////////////
// Read 2pt correlator ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// IDEA: First read in everything and compare checksums, if not correct give
// warning and exit
// afterwards look for correlator with tag. vectors tags and correlators have
// to be specified with the right size

void read_2pt_lime(const char* filename, std::vector<Tag>& tags,
                   std::vector<vec>& correlators){
  bool bigend = big_endian();
  if(tags.size() == correlators.size()){
    if(FILE* fp = fopen(filename, "r")){
      size_t global_check = 0;
      std::vector<size_t> checksums(tags.size());
      const int Lt = correlators[0].size();
      int MB_flag = 0; int ME_flag = 0;
      int file_status = 0;
      int act_status = 0;
      LimeReader* r = limeCreateReader(fp);
      n_uint64_t check_bytes = sizeof(size_t);
      n_uint64_t tag_bytes = sizeof(Tag);
      n_uint64_t data_bytes = correlators[0].size()*2*sizeof(double);
      file_status = limeReaderNextRecord(r);
      // From first message read global checksum
      MB_flag = limeReaderMBFlag(r);
      ME_flag = limeReaderMEFlag(r);
      //std::cout << MB_flag << " " << ME_flag << std::endl;
      if(MB_flag == 1 && ME_flag == 0){
        act_status = limeReaderReadData(&global_check, &check_bytes, r);
        if (bigend) global_check = swap_endian<size_t>(global_check);
        std::cout << global_check << std::endl;
      }
      //file_status = limeReaderNextRecord(r);
      // TODO: think about read in runinfo
      file_status = limeReaderNextRecord(r);
      // loop over all remaining messages
      int cnt = 0;
      do {
        //std::cout << "message: " << cnt/3 << std::endl;
        file_status = limeReaderNextRecord(r);
        MB_flag = limeReaderMBFlag(r);
        ME_flag = limeReaderMEFlag(r);
//        std::cout << MB_flag << " " << ME_flag << std::endl;
        // read correlator checksum
        if(MB_flag == 1 && ME_flag == 0){
          size_t check = 0;
          act_status = limeReaderReadData(&check, &check_bytes, r);
          if (bigend) check = swap_endian<size_t>(check);
          checksums.at(cnt/3) = check;
          //std::cout << check << std::endl;
        }
        // read correlator tag
        else if(MB_flag == 0 && ME_flag == 0){
          Tag read_tag;
          act_status = limeReaderReadData(&read_tag, &tag_bytes, r);
          if (bigend) read_tag = swap_single_tag(read_tag);
          tags.at(cnt/3) = read_tag;
        }
        // read correlator data
        else if(MB_flag == 0 && ME_flag == 1){
          vec corr;
          corr.resize(Lt);
          act_status = limeReaderReadData(corr.data(), &data_bytes, r);
          if (bigend) corr = swap_single_corr(corr);
          correlators.at(cnt/3) = corr;
        }
        ++cnt;
      } while (file_status != LIME_EOF);
      // Check file integrity
      size_t glob_bytes = correlators.size()*correlators[0].size();
      //concatenate all correlation functions in one vector
      std::vector<cmplx> collect;
      for(auto& c : correlators)
        for (auto& el : c) collect.push_back(el);
      file_check(global_check, checksums, collect);
    }
    else std::cout << "File does not exist!" << std::endl;
  }
  else std::cout << "#elements for tags and correlators not equal" << std::endl;
}

// Look for correlator corresponding to tag given as argument
void get_2pt_lime(const char* filename, const size_t num_corrs,
    const size_t corr_length, const Tag& tag,
    std::vector<cmplx >& corr){
  // Set up temporary structure as copy of incoming
  std::vector<vec> tmp_correlators(num_corrs);
  for(auto& corr : tmp_correlators) corr.resize(corr_length);
  std::vector<Tag> tmp_tags (num_corrs);
  // Read in whole Configuration function
  read_2pt_lime(filename, tmp_tags, tmp_correlators);
  // look for appropriate tag in vector of tags
  size_t tmp_tag_ind;
  for(size_t ind = 0; ind < tmp_tags.size(); ++ind)
    if (compare_tags(tag, tmp_tags[ind])) tmp_tag_ind = ind;
  
  // store correlation function in output data
  corr = tmp_correlators[tmp_tag_ind];

}

///////////////////////////////////////////////////////////////////////////////
// Dump to ASCII //////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Temporal extent has to be given
void ASCII_dump_corr(const char* filename, const size_t Lt,
                    const size_t num_corrs, const size_t p_sq,
                    const size_t l_sq) {
  // Set up temporary structure as copy of incoming
  std::vector<vec> correlators(num_corrs);
  for (auto& corr : correlators) corr.resize(Lt);
  std::vector<Tag> tags (num_corrs);

  // Read in whole Configuration function
  read_2pt_lime(filename, tags, correlators);
  std::vector<Tag> tags_found;
  std::vector<vec> corrs_found;

  // Append matching tags and correlators to vectors
  for(size_t n = 0; n < num_corrs; ++n){
    if (tags[n].mom_cm == p_sq) {
        std::array<int, 3> d1 = std_arr(tags[n].dis[1]);
        std::array<int, 3> d2 = std_arr(tags[n].dis[0]);
        if (square_comp(d1, d2) == l_sq){
        tags_found.push_back(tags[n]);
        corrs_found.push_back(correlators[n]);
      }
    }
  }
  size_t n = corrs_found.size();
  if( n == 0) std::cout << "No matching Correlators found!" <<
                              std::endl;
  else{
    std::cout << "# " << n << " correlators read in from file: " <<
              filename << std::endl;
    std::cout << "# cm_mom_sq dis_so_sq t real(C[t]) imag(C[t])" <<
              std::endl;
   for (size_t num = 0; num < n; ++num){
    for (size_t t = 0; t < Lt; ++t){
      std::cout << p_sq << " " << l_sq << " " << t << " " <<
                std::scientific << corrs_found[num][t].real() << " " <<
                corrs_found[num][t].imag() << std::endl;
    }
   }
  }

}

