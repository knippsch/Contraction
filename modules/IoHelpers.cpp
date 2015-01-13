#include "IoHelpers.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

///////////////////////////////////////////////////////////////////////////////
// Helper Routines ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// swap endianess of global data



//Tag handling/////////////////////////////////////////////////////////////////

// add two three momenta
std::array<int, 3 > add_mom(const std::array<int,3> p1 , const std::array<int,3> p2){
  std::array<int, 3 > sum;
  for(size_t i = 0; i < sum.size(); ++i){
    sum[i] = p1[i] + p2[i];
  }
return sum;

}

// Calculate p^2
int square_comp(const std::array<int, 3>& p1,
                const std::array<int, 3>& p2){
  int square = 0;
  for (size_t i = 0; i < 3; ++i){
    square += p1[i]*p2[i];  
  }
  return square;
}

// Compare two tags of correlation functions
bool compare_tags(const Tag& tag1, const Tag& tag2){
  bool flag = true;
  if (tag1.mom_cm != tag2.mom_cm) flag = false;
  if (memcmp(tag1.mom, tag2.mom, sizeof(tag1.mom)) != 0) flag = false;
  if (memcmp(tag1.dis, tag2.dis, sizeof(tag1.dis)) != 0) flag = false;
  if (memcmp(tag1.gam, tag2.gam, sizeof(tag1.gam)) != 0) flag = false;
  return flag;
}


// Set the tag from two operator structures
void set_tag(Tag& tag, const std::pair<size_t, size_t>& i){

  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();

  tag.mom[0] = square_comp(op_Corr[i.first].p3, op_Corr[i.first].p3);
  tag.mom[1] = tag.mom[0];
    for(size_t c = 0; c < 3; ++c){
      tag.dis[0][c] =op_Corr[i.first].dis3[c]; 
    }
    for(size_t c = 0; c < 3; ++c){
      tag.dis[1][c] =op_Corr[i.second].dis3[c]; 
    }

  tag.gam[0][0] = op_Corr[i.first].gamma[0];
  tag.gam[1][0] = op_Corr[i.second].gamma[0];

}

// Set the tag from two operator structures
void set_tag(Tag& tag, const std::array<size_t, 4>& i){

  const vec_pdg_Corr op_Corr = global_data->get_op_Corr();

  //TODO: use loops for that. I'm too tired
  std::array<int,3> cm_1 = add_mom(op_Corr[i[0]].p3, op_Corr[i[1]].p3);
  tag.mom_cm = square_comp(cm_1, cm_1);
  tag.mom[0] = square_comp(op_Corr[i[0]].p3, op_Corr[i[0]].p3);
  tag.mom[1] = square_comp(op_Corr[i[1]].p3, op_Corr[i[1]].p3);
  tag.mom[2] = square_comp(op_Corr[i[2]].p3, op_Corr[i[2]].p3);
  tag.mom[3] = square_comp(op_Corr[i[3]].p3, op_Corr[i[3]].p3);

  for(size_t q = 0; q < 4; ++q){
    for (size_t c = 0; c < 3; ++c){
      tag.dis[q][c] = op_Corr[i[q]].dis3[c];
    }
    for (size_t c = 0; c < 4; ++c){
      tag.gam[q][c] = op_Corr[i[q]].gamma[c];
    }
  }

}

// Map configurations of an infile to a tag
// void map_input_tag(const char* infile, Tag in){
//  
// // Get options from file to tag
// 
//   //one line of config file
//   std::string config_line;
//   //key and value to divide string to
//   std::string keyword;
//   std::string value;
//   
//   //three map class instances for int double and const char
//   std::map <std::string, int> integers;
//   std::map <std::string, double> floatings;
// 
//   //stream for config file
//   std::fstream config;
//   config.open(infile, std::fstream::in);
// 
//   //get each line of configuration into config_line
//   while (!std::getline(config, config_line).eof()){
//     if(config_line[0] != '#') {
//         //find first occurence of '=' in each line and divide strings
//         unsigned pos = config_line.find_first_of("=");
//         keyword = config_line.substr(0,(pos-1));
//         value = config_line.substr( (pos+2), std::string::npos);
// 
//         // String valued map variables
//         // Set momenta first
//         // Center of mass momenta momenta
//         if (keyword == "mom_sq_cm") {
//           in.mom_cm = (int)value;
//         }
//         if (keyword == "mom_1") {
//           in.mom[0] = (int)value;
//         }
//         if (keyword == "mom_2") {
//           in.mom[1] = (int)value;
//         }
//         if (keyword == "mom_3") {
//           in.mom[2] = (int)value;
//         }
//         if (keyword == "mom_4") {
//           in.mom[3] = (int)value;
//         }
//         // Gamma structures need splitting of the three values
//         if (keyword == "gamma_1" ) {
//           in.gam[0][0] = atoi(value.substr(0).c_str());
//           in.gam[0][1] = atoi(value.substr(2).c_str());
//           in.gam[0][2] = atoi(value.substr(4).c_str());
//           in.gam[0][3] = atoi(value.substr(6).c_str());
//         }
//         if (keyword == "gamma_2" ) {
//           in.gam[1][0] = atoi(value.substr(0).c_str());
//           in.gam[1][1] = atoi(value.substr(2).c_str());
//           in.gam[1][2] = atoi(value.substr(4).c_str());
//           in.gam[1][3] = atoi(value.substr(6).c_str());
//         }
//         if (keyword == "gamma_3" ) {
//           in.gam[2][0] = atoi(value.substr(0).c_str());
//           in.gam[2][1] = atoi(value.substr(2).c_str());
//           in.gam[2][2] = atoi(value.substr(4).c_str());
//           in.gam[2][3] = atoi(value.substr(6).c_str());
//         }
//         if (keyword == "gamma_4") {
//           in.gam[3][0] = atoi(value.substr(0).c_str());
//           in.gam[3][1] = atoi(value.substr(2).c_str());
//           in.gam[3][2] = atoi(value.substr(4).c_str());
//           in.gam[3][3] = atoi(value.substr(6).c_str());
//         }
// 
//         // Displacement vectors
//         if (keyword == "dis_1"){
//           in.dis[0][0] = ;
//           in.dis[0][1] = ;
//           in.dis[0][2] = ;
//         }
//         if (keyword == "dis_2"){
//           in.dis[1][0] = ;
//           in.dis[1][1] = ;
//           in.dis[1][2] = ;
//         }
//         if (keyword == "dis_3"){
//           in.dis[2][0] = ;
//           in.dis[2][1] = ;
//           in.dis[2][2] = ;
//         }
//         if (keyword == "dis_4"){
//           in.dis[3][0] = ;
//           in.dis[3][1] = ;
//           in.dis[3][2] = ;
//         }
//     }
//   }
//   config.close();
// 
// // CM momentum squared
//  in.mom_cm = ;
// 
// // Source
// // 1st particle
//  in.mom[0] = ;
//  in.dis[0][0] = ;
//  in.dis[0][1] = ;
//  in.dis[0][2] = ;
// 
//  in.gam[0][0] = ;
//  in.gam[0][1] = ;
//  in.gam[0][2] = ;
//  in.gam[0][3] = ;
// 
//  // 2nd particle 
//  in.mom[1] = ;
//  in.dis[1][0] = ;
//  in.dis[1][1] = ;
//  in.dis[1][2] = ;
// 
//  in.gam[1][0] = ;
//  in.gam[1][1] = ;
//  in.gam[1][2] = ;
//  in.gam[1][3] = ;
// 
//  // Sink
//  // 1st particle
//  in.mom[2] = ;
//  in.dis[2][0] = ;
//  in.dis[2][1] = ;
//  in.dis[2][2] = ;
// 
//  in.gam[2][0] = ;
//  in.gam[2][1] = ;
//  in.gam[2][2] = ;
//  in.gam[2][3] = ;
// 
//  // 2nd particle
//  in.mom[3] = ;
//  in.dis[3][0] = ;
//  in.dis[3][1] = ;
//  in.dis[3][2] = ;
// 
//  in.gam[3][0] = ;
//  in.gam[3][1] = ;
//  in.gam[3][2] = ;
//  in.gam[3][3] = ;
// }

//TODO
// Convert ascii labels to correlation tag
//Tag id(size_t g_so, size_t g_si, size_t p_so, size_t p_si, size_t dis_so, size_t dis_si){
//  Tag item;
//  item.mom[0] = p_so;
//  item.mom[1] = p_si;
//  item.dis[0] = dis_so;
//  item.dis[1] = dis_si;
//  item.gam[0] = g_so;
//  item.gam[1] = g_si;
//  return item;
//} 

// File handling ///////////////////////////////////////////////////////////////

// write global checksum, and runinfo in new file
void write_1st_msg(const char* filename, GlobalDat& dat,
                          size_t chksum){
  // open file
  FILE* fp;
  fp = fopen(filename, "w");
  LimeWriter* w;
  w = limeCreateWriter(fp);
  // first record is global checksum
  int MB_flag = 1;
  int ME_flag = 0;
  LimeRecordHeader* chk;
  n_uint64_t chk_bytes = sizeof(chksum);
  chk = limeCreateHeader(MB_flag, ME_flag, "Global Checksum", chk_bytes);
  limeWriteRecordHeader(chk, w);
  limeDestroyHeader(chk);
  limeWriteRecordData(&chksum, &chk_bytes, w);
  
  // second record is data on calculation
  MB_flag = 0; ME_flag = 1;
  
  LimeRecordHeader* run_id;
  n_uint64_t run_id_bytes = sizeof(dat);
  run_id = limeCreateHeader(MB_flag, ME_flag, "Runinfo", run_id_bytes);
  limeWriteRecordHeader(run_id,w);
  limeDestroyHeader(run_id);
  limeWriteRecordData(&dat, &run_id_bytes, w);
  limeDestroyWriter(w);
  fclose(fp);
} 

static bool tag_exist(const char* filename, Tag id){
  bool found = false;
  
  return found;
}

//append correlation functions to file
void append_msgs(const char* filename, std::vector<vec>& corr, std::vector<Tag>& tags,
              LimeWriter* w, bool be){
  // Each message contains three records:
  // 1st: Checksum for the Correlator
  // 2nd: Tag for the Correlator
  // 3rd: Correlationfunction itself
  LimeRecordHeader* corr_chk;
  LimeRecordHeader* id;
  LimeRecordHeader* corr_hd;

  boost::uint64_t corr_chksum;
  n_uint64_t tag_bytes = sizeof(tags[0]);
  n_uint64_t data_bytes = (corr[0]).size()*2*sizeof(double);
  char headername[100];  
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  for(size_t el = 0; el < corr.size(); ++el){

    // 1st record for checksum
    corr_chksum = checksum <vec> (corr[el], corr[el].size());
    if(be) corr_chksum = swap_endian<boost::uint64_t>(corr_chksum);
    n_uint64_t chk_bytes = sizeof(corr_chksum);
    ME_flag = 0; MB_flag = 1;
    corr_chk = limeCreateHeader(MB_flag, ME_flag,"Correlator checksum", chk_bytes);
    limeWriteRecordHeader( corr_chk, w );
    limeDestroyHeader( corr_chk );
    limeWriteRecordData( &corr_chksum, &chk_bytes, w );

    // 2nd record for Tag as struct of 3 2dim integer-arrays
    ME_flag = 0; MB_flag = 0;
    sprintf(headername,"Tag of Correlator with p^2 = %zd", el);
    id = limeCreateHeader( MB_flag, ME_flag, headername , tag_bytes );
    limeWriteRecordHeader( id, w );
    limeDestroyHeader( id );
    limeWriteRecordData( &tags[el], &tag_bytes, w );

    // 3rd record for correlator belonging to tag
    ME_flag = 1; MB_flag = 0; 
    sprintf(headername,"Correlator with p^2 = %zd", el);
    corr_hd = limeCreateHeader( MB_flag, ME_flag, headername, data_bytes );
    limeWriteRecordHeader( corr_hd, w );
    limeDestroyHeader( corr_hd );
    limeWriteRecordData( corr[el].data(), &data_bytes, w ); 
  }
}

// Check checksums
void file_check(const size_t glob_check,
                const std::vector<boost::uint64_t>& checksums,
                const std::vector<cmplx>& correlators){
  size_t tmp = 0;
  size_t bytecount = correlators.size()*2*sizeof(double);
  
  tmp = checksum<std::vector<cmplx> > (correlators, bytecount);
  if (big_endian()) swap_endian <size_t> (tmp);
  std::cout << bytecount << std::endl;
  if (tmp == glob_check) {
    std::cout << " File Checksum matches data. "<< std::endl;
  }
  else{
    std::cout << " Checksum broken, please check correlation functions ! "
       << std::endl;
    std::cout << "Read in: " << glob_check << " calculated: " << tmp << std::endl; 
  }

}


