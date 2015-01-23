#include "CorrelatorIo2pt.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//TODO: put that into the IoHelpers. Does not work naively as the template
// type must be known at compile time

// convert multiarray 2pt correlator to vector to match write_2pt_lime

 //TODO: functions need to be adjusted to new tag structure 
// Set the tag from two operator structures


void set_tag(Tag& tag, const std::pair<size_t, size_t>& i, const vec_index_2pt& corr_type){

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  
  tag.q_cont = "udud";

  tag.mom.push_back(op_Corr[corr_type[i.first].index_Q2].p3);
  tag.mom.push_back(op_Corr[corr_type[i.first].index_Corr].p3);
  tag.mom.push_back(op_Corr[corr_type[i.second].index_Q2].p3);
  tag.mom.push_back(op_Corr[corr_type[i.second].index_Corr].p3);
  
  tag.mom_cm = square_comp(add_mom(tag.mom[0], tag.mom[2]),
                           add_mom(tag.mom[0], tag.mom[2]));

  tag.dis.push_back(op_Corr[corr_type[i.first].index_Q2].dis3);
  tag.dis.push_back(op_Corr[corr_type[i.first].index_Corr].dis3);
  tag.dis.push_back(op_Corr[corr_type[i.second].index_Q2].dis3);
  tag.dis.push_back(op_Corr[corr_type[i.second].index_Corr].dis3);

  std::array<int, 4> tmp = {{0,0,0,0}};
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i.first].index_Q2].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i.first].index_Q2].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i.first].index_Corr].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i.first].index_Corr].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i.first].index_Q2].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i.second].index_Q2].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i.first].index_Corr].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i.second].index_Corr].gamma[ind];
  }
  tag.gam.push_back(tmp);

}

// Set the tag from two operator structures
void set_tag(Tag& tag, const size_t i, const vec_index_4pt& corr_type){

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  tag.q_cont = "udud";

  tag.mom.push_back(op_Corr[corr_type[i].index_Q2[0]].p3);
  tag.mom.push_back(op_Corr[corr_type[i].index_Corr[0]].p3);
  tag.mom.push_back(op_Corr[corr_type[i].index_Q2[1]].p3);
  tag.mom.push_back(op_Corr[corr_type[i].index_Corr[1]].p3);
  
  tag.mom_cm = square_comp(add_mom(tag.mom[0], tag.mom[2]),
                           add_mom(tag.mom[0], tag.mom[2]));

  tag.dis.push_back(op_Corr[corr_type[i].index_Q2[0]].dis3);
  tag.dis.push_back(op_Corr[corr_type[i].index_Corr[0]].dis3);
  tag.dis.push_back(op_Corr[corr_type[i].index_Q2[1]].dis3);
  tag.dis.push_back(op_Corr[corr_type[i].index_Corr[1]].dis3);

  std::array<int, 4> tmp = {{0,0,0,0}};
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Q2[0]].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Q2[0]].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Corr[0]].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Corr[0]].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Q2[1]].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Q2[1]].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Corr[1]].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Corr[1]].gamma[ind];
  }
  tag.gam.push_back(tmp);

}

// Set the tag from two operator structures
void set_tag(Tag& tag, const size_t i, const vec_index_2pt& corr_type){

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  tag.q_cont = "ud";

  tag.mom.push_back(op_Corr[corr_type[i].index_Q2].p3);
  tag.mom.push_back(op_Corr[corr_type[i].index_Corr].p3);
  
  tag.mom_cm = square_comp(tag.mom[0], tag.mom[0]);

  tag.dis.push_back(op_Corr[corr_type[i].index_Q2].dis3);
  tag.dis.push_back(op_Corr[corr_type[i].index_Corr].dis3);

  std::array<int, 4> tmp = {{0,0,0,0}};
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Q2].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Q2].gamma[ind];
  }
  tag.gam.push_back(tmp);
  for(auto& el : tmp) el = 0;
  for(size_t ind = 0; 
      ind < op_Corr[corr_type[i].index_Corr].gamma.size(); ind++){
    tmp[ind] = op_Corr[corr_type[i].index_Corr].gamma[ind];
  }
  tag.gam.push_back(tmp);

}

void set_tag(Tag& tag, const std::pair<size_t, size_t>& i, 
             const vec_index_4pt& corr_type){

  //This should not be called
  exit(0);
}

template <typename io_list> 
void convert_hadron_to_vec(const io_list& op_io, const array_cd_d2& C2_mes,
                           const std::string& corr_type, std::vector<Tag>& tags,
                           std::vector<vec>& corr){

  const size_t Lt = global_data->get_Lt();

  vec_index_2pt lookup_2pt = global_data->get_lookup_2pt_trace(); 
  vec_index_4pt lookup_4pt = global_data->get_lookup_4pt_trace();

  corr.resize(op_io.size());
  for(auto& c : corr)
    c.resize(Lt);
  tags.resize(op_io.size());

  for(const auto& op : op_io){
  for(const auto& i : op.index_pt){
    //TODO: solve the copying of C2_mes (boost) into corr (std::vec) without
    // the copy constructor in assign()
    corr[op.id].assign(C2_mes[op.id].origin(), C2_mes[op.id].origin() + 
        C2_mes[op.id].size());

    if( (corr_type == "C2+") || (corr_type == "C4I2+_1") || (corr_type == "C4I2+_2") ){
      set_tag(tags[op.id], i, lookup_2pt);
    }
    else if( corr_type == "C4I2+_3"){
      set_tag(tags[op.id], i, lookup_4pt);
    }
  }}
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void export_corr_IO (const char* filename, const vec_index_IO_1& op_IO,
                     const std::string& corr_type, const array_cd_d2& C2_mes){
  GlobalDat dat;
  std::vector<Tag> tags;
  std::vector<std::string> tag_strings;
  std::string tag_string;
  std::vector<vec> corr;

  convert_hadron_to_vec <vec_index_IO_1> (op_IO, C2_mes, corr_type, tags, corr);
  for(const auto& tag : tags){
    tag_to_string(tag, tag_string);
    tag_strings.push_back(tag_string);
  }

//  if (file_exist(filename)){
//    char filename_new [150];
//    sprintf(filename_new,"_changed");
//    std::cout << "file already exists! Renaming to "
//      << filename_new << std::endl;
//  }
  std::cout << "\twriting " << corr_type << " to " << filename << std::endl;
  write_2pt_lime(filename, dat, tag_strings, corr);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void export_corr_IO (const char* filename, const vec_index_IO_2& op_IO,
                     const std::string& corr_type, const array_cd_d2& C2_mes){
  GlobalDat dat;
  std::vector<Tag> tags;
  std::vector<std::string> tag_strings;
  std::string tag_string;
  std::vector<vec> corr;

  convert_hadron_to_vec <vec_index_IO_2> (op_IO, C2_mes,corr_type, tags, corr);
  for(const auto& tag : tags){
    tag_to_string(tag, tag_string);
    tag_strings.push_back(tag_string);
  }

//  if (file_exist(filename)){
//    char filename_new [150];
//    sprintf(filename_new,"_changed");
//    std::cout << "file already exists! Renaming to "
//      << filename_new << std::endl;
//  }
  std::cout << "\twriting " << corr_type << " to " << filename << std::endl;
  write_2pt_lime(filename, dat, tag_strings, corr);
}

///////////////////////////////////////////////////////////////////////////////
// Write 2pt correlator ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void write_2pt_lime(const char* filename, GlobalDat& dat, 
    std::vector<std::string>& tags, std::vector<vec>& corr){

  bool be = big_endian();
  if (be){
    swap_correlators(corr);
    //swap_tag_vector(tags);
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
// Read npt correlator ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// IDEA: First read in everything and compare checksums, if not correct give
// warning and exit
// afterwards look for correlator with tag. vectors tags and correlators have
// to be specified with the right size

void read_2pt_lime(const char* filename, std::vector<std::string>& tags,
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
      n_uint64_t tag_bytes = tags[0].length()*sizeof(char);
      std::cout << tag_bytes << std::endl;
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
          // std::cout << check << std::endl;
        }
        // read correlator tag
        else if(MB_flag == 0 && ME_flag == 0){
          std::string read_tag;
          read_tag.resize(tags[0].length());
          act_status = limeReaderReadData(&read_tag[0], &tag_bytes, r);
          //std::cout << read_tag << std::endl;
          //if (bigend) read_tag = swap_single_tag(read_tag);
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
    const size_t corr_length, const std::string& tag,
    std::vector<cmplx >& corr){
  // Set up temporary structure as copy of incoming
  std::vector<vec> tmp_correlators(num_corrs);
  for(auto& corr : tmp_correlators) corr.resize(corr_length);
  std::vector<std::string> tmp_tags (num_corrs);
  // Read in whole Configuration function
  read_2pt_lime(filename, tmp_tags, tmp_correlators);
  // look for appropriate tag in vector of tags
  size_t tmp_tag_ind;
  for(size_t ind = 0; ind < tmp_tags.size(); ++ind)
    if (tag == tmp_tags[ind]) tmp_tag_ind = ind;

  // store correlation function in output data
  corr = tmp_correlators[tmp_tag_ind];

}

///////////////////////////////////////////////////////////////////////////////
// Dump to ASCII //////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Temporal extent has to be given
void ASCII_dump_corr(const char* filename, const char* infile, const size_t Lt,
    const size_t num_corrs) {
  // Set up temporary structure as copy of incoming
  std::vector<vec> correlators(num_corrs);
  for (auto& corr : correlators) corr.resize(Lt);
  std::vector<std::string> tags (num_corrs);
  std::string query;
  std::string query_one;
  Tag id;
  compose_string(infile, query);
  string_to_tag(query, id);
  tag_to_string(id, query_one);
  std::cout << "Looking for: " << std::endl;
  std::cout << query_one << std::endl;
  for(auto& el : tags) el.resize(query_one.length());
  // Read in whole Configuration function
  read_2pt_lime(filename, tags, correlators);
  std::cout << "finished read in" << std::endl;
  std::vector<std::string> tags_found;
  std::vector<vec> corrs_found;

  // Append matching tags and correlators to vectors
  for(size_t n = 0; n < num_corrs; ++n){
    if(tags[n] == query_one){
      tags_found.push_back(tags[n]);
      corrs_found.push_back(correlators[n]);
    }
  }
  size_t n = corrs_found.size();
  if( n == 0) std::cout << "No matching Correlators found!" <<
    std::endl;
  else{
    std::cout << "# " << n << " correlators read in from file: " <<
      filename << std::endl;
    std::cout << "# t real(C[t]) imag(C[t])" <<
      std::endl;
    for (size_t num = 0; num < n; ++num){
      for (size_t t = 0; t < Lt; ++t){
        std::cout << t << " " <<
          std::scientific << corrs_found[num][t].real() << " " <<
          corrs_found[num][t].imag() << std::endl;
      }
    }
  }

}

