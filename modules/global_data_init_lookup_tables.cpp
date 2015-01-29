#include "global_data.h"
#include "global_data_utils.h"

namespace gdu = ::global_data_utils;

namespace {

// need some functions from namespace global_data_utils
using gdu::add_p3;
using gdu::abs_p3;
using gdu::compare_quantum_numbers_of_pdg;
using gdu::compare_mom_dis_of_pdg;
using gdu::compare_index_list;
using gdu::set_index_corr;
using gdu::set_index_2pt;
using gdu::set_index_4pt;

void init_lookup_corr(const Correlator_list& correlator_list, 
                      const std::vector<Operator_list>& operator_list,
                      vec_pdg_Corr& lookup_corr, vec_pd_VdaggerV& lookup_vdv,
                      vec_pd_rVdaggerVr& lookup_rvdvr, size_t index_of_unity) {

  // extracting all operators which are used in correlations functions
  std::vector<int> used_operators;
  for(const auto& corr_list : correlator_list)
    used_operators.insert(used_operators.end(), 
                          corr_list.operator_numbers.begin(), 
                          corr_list.operator_numbers.end());
  
  sort(used_operators.begin(), used_operators.end());
  used_operators.erase(std::unique(used_operators.begin(), 
                                   used_operators.end()),
                       used_operators.end());
  // write quantum number in lookup_corr
  for(const auto& op_entry : used_operators){
    for(const auto& individual_operator : operator_list[op_entry]){
      pdg write;
      write.gamma = individual_operator.gammas;
      write.dis3 = individual_operator.dil_vec;
      for(const auto& mom_vec : individual_operator.mom_vec){
        for(auto mom : mom_vec){
          write.p3 = mom;
          lookup_corr.push_back(write);
        }
      }
    }
  }
  // doubly counted lookup_corr entries are deleted
  auto it = lookup_corr.begin();
  while(it != lookup_corr.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_corr.end()) {
      if(compare_quantum_numbers_of_pdg(*it, *it2))
        lookup_corr.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // sorting lookup_corr for equal momentum and displacement vectors - makes it
  // easier to run over it with auto loops
  std::vector<pdg> dump_write;
  while(lookup_corr.size() != 0){

    it = lookup_corr.begin();
    dump_write.push_back(*it);
    lookup_corr.erase(it);

    auto it2 = dump_write.end()-1;
    while(it != lookup_corr.end()) {
      if(compare_mom_dis_of_pdg(*it, *it2)){
        dump_write.push_back(*it);
        lookup_corr.erase(it);
      }
      else
        it++;
    }
  }
  lookup_corr.swap(dump_write);


  // setting the identification numbers of lookup_corr
  size_t counter = 0;
  for(auto& op : lookup_corr){
    op.id = counter++;
  }

  // final setting lookuptables for vdaggerv and so on
  set_index_corr(lookup_corr, lookup_vdv, lookup_rvdvr);

  // setting index_of_unity
  pdg zero;
  for(auto& mom : zero.p3)
    mom = 0;
  for(auto& dis : zero.dis3)
    dis = 0;
  for(const auto& op : lookup_corr)
    if( (op.first_vdv == true) && compare_mom_dis_of_pdg(zero, op)){
      index_of_unity = op.id;
      break;
    }

  std::cout << "lookup_corr" << std::endl;
  for(auto a : lookup_corr){
    std::cout << a.id;
    for(auto b : a.gamma)
    std::cout << " " << b;
    for(auto b : a.dis3)
      std::cout << " " << b;
    for(auto b : a.p3)
      std::cout << " " << b;
    std::cout << "\t" << a.id_vdv << " " << a.first_vdv << " " 
              << a.negative_momentum << " " << a.id_rvdvr << std::endl;
  }

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_2pt(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     const vec_pdg_Corr& lookup_corr, 
                     vec_index_2pt& lookup_2pt) {

  // TODO: think about symmetry of switching index_Q2 and index_Corr
  // init lookup_2pt
  for(const auto& corr : correlator_list){

    // must be done for 2pt as well as 4pt function
    if( (corr.type.compare(0,3,"C2+") == 0) || 
        (corr.type.compare(0,5,"C4I2+") == 0) ){
      for(const auto& infile_op_so : operator_list[corr.operator_numbers[0]]){
      for(const auto& infile_op_si : operator_list[corr.operator_numbers[1]]){
        // TODO: give that guy the quark numbers from corr and think about 
        // giving the random numbers as well
        set_index_2pt(infile_op_so, infile_op_si, lookup_corr, lookup_2pt);
      }} //loops over same physical situation end here
    } //case 2pt-function ends here

    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& infile_op_so : operator_list[corr.operator_numbers[2]]){
      for(const auto& infile_op_si : operator_list[corr.operator_numbers[3]]){
        set_index_2pt(infile_op_so, infile_op_si, lookup_corr, lookup_2pt);
      }} //loops over same physical situation end here
    } //case 2pt-function ends here

  } // loop over correlator_list ends here

  // doubly counted lookup_index_C2 entries are deleted in order to not be
  // calculated twice
  auto it = lookup_2pt.begin();
  while(it != lookup_2pt.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_2pt.end()) {
      if( (it->index_Q2 == it2->index_Q2) && 
          (it->index_Corr == it2->index_Corr) )
        lookup_2pt.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // initialize the id's of lookup_2pt
  size_t counter = 0;
  for(auto& op : lookup_2pt){
    op.id = counter++;
  }

  std::cout << "lookup_2pt" << std::endl;
  for(const auto& a : lookup_2pt){
    std::cout << a.id << "\t" << a.index_Q2 << "\t" << a.index_Corr<< std::endl;
  }
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_C2plus_IO(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     const vec_pdg_Corr& lookup_corr, vec_index_2pt& lookup_2pt,
                     vec_index_IO_1& lookup_2pt_IO) {

  std::array<int, 3> zero = {{0, 0, 0}};

  // init lookup_2pt_IO
  for(const auto& corr : correlator_list){
    if(corr.type.compare(0,3,"C2+") == 0){
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
        for(auto op : lookup_2pt){
          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2], 
                                            op1_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr], 
                                            op2_from_list)){
            // momentum conservation in case of 2pt function including momentum turning
            // due to gamma_5 trick
            if(add_p3(lookup_corr[op.index_Q2], lookup_corr[op.index_Corr]) == zero){
              index_IO_1 write;
              // write only one value in each indexlist
              write.index_pt.emplace_back(op.id);
              lookup_2pt_IO.push_back(write);
            }
          }}
        }
      }}
    }
  }

  // doubly counted lookup_corr entries are deleted
  auto it = lookup_2pt_IO.begin();
  while(it != lookup_2pt_IO.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_2pt_IO.end()) {
      if(compare_index_list(*it, *it2))
        lookup_2pt_IO.erase(it2);
      else
        it2++;
    }
    it++;
  }

  size_t counter = 0;
  for(auto& op_IO : lookup_2pt_IO){
    op_IO.id = counter;
    counter++;
  }

  std::cout << "lookup_2pt_IO" << std::endl;
  for(const auto& a : lookup_2pt_IO){
    for(const auto& b : a.index_pt)
      std::cout << a.id << "\t" << b << std::endl;
    std::cout << std::endl;
  }
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_C4I2plus_IO(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     vec_index_2pt& lookup_2pt, vec_index_4pt& lookup_4pt,
                     vec_index_IO_2& lookup_4pt_1_IO, 
                     vec_index_IO_2& lookup_4pt_2_IO) {

  std::array<int, 3> zero = {{0,0,0}};

  // init lookup_4pt_IO
  for(const auto& corr : correlator_list){
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& op_C4 : lookup_4pt){
        for(const auto& op : lookup_2pt){
        if( (op_C4.index_Q2[0] == op.index_Q2) && (op_C4.index_Corr[0] == op.index_Corr)){
          for(const auto& op2 : lookup_2pt){
          if( (op_C4.index_Q2[1] == op2.index_Q2) && (op_C4.index_Corr[1] == op2.index_Corr)){
            index_IO_2 write;
            write.index_pt.emplace_back
                (std::pair<size_t, size_t>(op.id, op2.id));
            lookup_4pt_1_IO.push_back(write);
            lookup_4pt_2_IO.push_back(write);
          }}
        }}
      }
    }
  }

//      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
//      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
//      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
//      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
//        for(auto op : lookup_2pt){
//          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
//          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2], 
//                                            op1_from_list)){
//          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr], 
//                                            op2_from_list)){
//          for(auto op2 : lookup_2pt){
//            if(compare_quantum_numbers_of_pdg(lookup_corr[op2.index_Q2], 
//                                              op3_from_list)){
//            if(compare_quantum_numbers_of_pdg(lookup_corr[op2.index_Corr], 
//                                              op4_from_list)){
//              
//              if(add_p3(lookup_corr[op.index_Q2], lookup_corr[op2.index_Q2]) 
//                    == zero){
//              if(abs_p3(lookup_corr[op.index_Q2]) == abs_p3(lookup_corr[op2.index_Q2])){
//                index_IO_2 write;
//                write.index_pt.emplace_back
//                    (std::pair<size_t, size_t>(op.id, op2.id));
//                lookup_4pt_1_IO.push_back(write);
//                lookup_4pt_2_IO.push_back(write);
//              }}
//            }}
//          }
//          }}
//        }
//      }}}}
//    }
//  }

  // doubly counted lookup_corr entries are deleted
  auto it1 = lookup_4pt_1_IO.begin();
  while(it1 != lookup_4pt_1_IO.end()) {
    auto it2 = it1;
    it2++;
    while(it2 != lookup_4pt_1_IO.end()) {
      if(compare_index_list(*it1, *it2))
        lookup_4pt_1_IO.erase(it2);
      else
        it2++;
    }
    it1++;
  }

  size_t counter = 0;
  for(auto& op_IO : lookup_4pt_1_IO){
    op_IO.id = counter;
    counter++;
  }

  // doubly counted lookup_corr entries are deleted
  auto it3 = lookup_4pt_2_IO.begin();
  while(it3 != lookup_4pt_2_IO.end()) {
    auto it4 = it3;
    it4++;
    while(it4 != lookup_4pt_2_IO.end()) {
      if(compare_index_list(*it3, *it4))
        lookup_4pt_2_IO.erase(it4);
      else
        it4++;
    }
    it3++;
  }

  counter = 0;
  for(auto& op_IO : lookup_4pt_2_IO){
    op_IO.id = counter;
    counter++;
  }

  std::cout << "lookup_4pt_1_IO" << std::endl;
  for(const auto& a : lookup_4pt_1_IO){
    for(const auto& b : a.index_pt)
      std::cout << a.id << "\t" << b.first << "\t" << b.second << std::endl;
    std::cout << std::endl;
  }

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_4pt(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     const vec_pdg_Corr& lookup_corr, vec_index_4pt& lookup_4pt,
                     vec_index_IO_1& lookup_4pt_3_IO) {

  // initialization of op_C4
  
  for(auto& corr : correlator_list){
    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        set_index_4pt(op1_from_list, op2_from_list, op3_from_list, 
                      op4_from_list, lookup_corr, lookup_4pt);
      }}}} //loops over same physical situation end here
    } //case 4pt-function ends here
  } // loop over correlator_list ends here

  // doubly counted lookup_index_C4 entries are deleted in order to not be
  // calculated twice
  auto it = lookup_4pt.begin();
  while(it != lookup_4pt.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_4pt.end()) {
      if( (it->index_Q2[0] == it2->index_Q2[0]) && 
          (it->index_Corr[0] == it2->index_Corr[0]) &&
          (it->index_Q2[1] == it2->index_Q2[1]) && 
          (it->index_Corr[1] == it2->index_Corr[1]) )
        lookup_4pt.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // initialize the id's of lookup_2pt
  size_t counter = 0;
  for(auto& op : lookup_4pt){
    op.id = counter++;
  }

  std::cout << "lookup_4pt" << std::endl;
  for(auto a : lookup_4pt){
    std::cout << a.id << "\t" << a.index_Q2[0] << "\t" << a.index_Corr[0] 
              << "\t" << a.index_Q2[1] << "\t" << a.index_Corr[1] << std::endl;
  }

  // init lookup_4pt_IO
  for(const auto& corr : correlator_list){

    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        for(auto op : lookup_4pt){
          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2[0]], 
                                            op1_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr[0]], 
                                            op2_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2[1]], 
                                            op3_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr[1]], 
                                            op4_from_list)){
            index_IO_1 write;
            write.index_pt.emplace_back(op.id);
            lookup_4pt_3_IO.push_back(write);
          }}}}
        }
      }}}}
    }
  }

  counter = 0;
  for(auto& op : lookup_4pt_3_IO){
    op.id = counter++;
  }

  std::cout << "lookup_4pt_3_IO" << std::endl;
  for(const auto& a : lookup_4pt_3_IO){
    for(const auto& b : a.index_pt)
      std::cout << "\t" << a.id << "\t" << b << std::endl;
    std::cout << std::endl;
  }

}

// function to obtain the index combinations of the random vectors
// for one quark
void set_rnd_vec_1pt(std::vector<quark>& quarks, indexlist_1& rnd_vec_1pt) {
  // ATM hardcoded, check which quarks are being used
  const int q1 = 0;
  const int rndq1 = quarks[q1].number_of_rnd_vec;

  // check if there are enough random vectors
  if(rndq1 < 1) {
    std::cerr << "There are not enough random vectors for 1point functions" << std::endl;
    exit(-1);
  }
  for(size_t i = 0; i < rndq1; ++i) {
    rnd_vec_1pt.emplace_back(i);
  }

//  std::cout << "rnd_vec test: " << rnd_vec_C1.size() << std::endl;
//  for(auto& r : rnd_vec_C1) {
//    std::cout << r << std::endl;
//  }
}

// function to obtain the index combinations of the random vectors
// for two quarks
void set_rnd_vec_2pt(std::vector<quark>& quarks, indexlist_2& rnd_vec_2pt) {
  // ATM hardcoded, check which quarks are being used
  const int q1 = 0, q2 = 0;
  const int rndq1 = quarks[q1].number_of_rnd_vec;
  const int rndq2 = quarks[q2].number_of_rnd_vec;

  // check if there are enough random vectors
  if( (q1 == q2) && (rndq1 < 2) ) {
    std::cerr << "There are not enough random vectors for 2point functions" << std::endl;
    exit(-1);
  }
  for(size_t i = 0; i < rndq1; ++i) {
    for(size_t j = 0; j < rndq1; ++j) {
      // if the quarks are the same, skip the combinations
      // with i==j
      if( (q1 != q2) || (i != j) ) {
        rnd_vec_2pt.emplace_back(i, j);
      }
    }
  }

//  std::cout << "rnd_vec test: " << rnd_vec_C2.size() << std::endl;
//  for(auto& r : rnd_vec_C2) {
//    std::cout << r.first << " " << r.second << std::endl;
//  }
}

// function to obtain index combinations of the random vectors 
// for three quarks
void set_rnd_vec_3pt(std::vector<quark>& quarks, indexlist_3& rnd_vec_3pt) {
  // ATM hardcoded, check which quarks are being used
  const int q1 = 0;
  const int q2 = 0;
  const int q3 = 0;
  const int rndq1 = quarks[q1].number_of_rnd_vec;
  const int rndq2 = quarks[q2].number_of_rnd_vec;
  const int rndq3 = quarks[q3].number_of_rnd_vec;
  //check if there are enough random vectors
  if( (q1 == q2 ) && (rndq1 < 2) ) {
    if( ( (q1 == q3) && (q2 == q3) && (rndq1 < 3) ) ||
        ( (q1 == q3) && (q2 != q3) && (rndq1 < 2) ) ||
        ( (q1 != q3) && (q2 == q3) && (rndq2 < 2) ) ) {
      std::cerr << "there are not enough random vectors for 4point functions" << std::endl;
      exit(-1);
    }
  }

  // Build all combinations of random vectors for 4point functions.
  // if two or more quarks have the same id, skip combinations where at least
  // two are equal
  for(size_t rnd1 = 0; rnd1 < rndq1; ++rnd1) {
    for(size_t rnd2 = 0; rnd2 < rndq2; ++rnd2) {
      if( (q1 != q2 ) || ( (q1 == q2) && (rnd1 != rnd2) ) ) {
        for(size_t rnd3 = 0; rnd3 < rndq3; ++rnd3) {
          if( ( (q1 != q3) && (q2 != q3) ) ||
              ( (q1 == q3) && (q2 == q3) && (rnd1 != rnd3) && (rnd2 != rnd3) ) ||
              ( (q1 == q3) && (q2 != q3) && (rnd1 != rnd3) ) ||
              ( (q1 != q3) && (q2 == q3) && (rnd2 != rnd3) ) ) {
            rnd_vec_3pt.emplace_back(std::array<size_t, 3> {{rnd1, rnd2, rnd3}});
          }
        }
      }
    }
  }

//  std::cout << "rnd_vec test: " << rnd_vec_C3.size() << std::endl;
//  for(auto& r : rnd_vec_C3) {
//    std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
//  }
}

// function to obtain index combinations of the randomvectors for four quarks
void set_rnd_vec_4pt(std::vector<quark>& quarks, indexlist_4& rnd_vec_4pt) {
  // ATM hardcoded, check which quarks are being used
  const int q1 = 0;
  const int q2 = 0;
  const int q3 = 0;
  const int q4 = 0;
  const int rndq1 = quarks[q1].number_of_rnd_vec;
  const int rndq2 = quarks[q2].number_of_rnd_vec;
  const int rndq3 = quarks[q3].number_of_rnd_vec;
  const int rndq4 = quarks[q4].number_of_rnd_vec;
  //check if there are enough random vectors
  if( (q1 == q2 ) && (rndq1 < 2) ) {
    if( ( (q1 == q3) && (q2 == q3) && (rndq1 < 3) ) ||
        ( (q1 == q3) && (q2 != q3) && (rndq1 < 2) ) ||
        ( (q1 != q3) && (q2 == q3) && (rndq2 < 2) ) ) {
      if( ( (q1 == q4) && (q2 == q4) && (q3 == q4) && (rndq1 < 4) ) ||
          ( (q1 == q4) && (q2 != q4) && (q3 != q4) && (rndq1 < 2) ) ||
          ( (q1 != q4) && (q2 == q4) && (q3 != q4) && (rndq2 < 2) ) ||
          ( (q1 != q4) && (q2 != q4) && (q3 == q4) && (rndq3 < 2) ) ||
          ( (q1 == q4) && (q2 == q4) && (q3 != q4) && (rndq1 < 3) ) ||
          ( (q1 == q4) && (q2 != q4) && (q3 == q4) && (rndq1 < 3) ) ||
          ( (q1 != q4) && (q2 == q4) && (q3 == q4) && (rndq2 < 3) ) ) {
        std::cerr << "there are not enough random vectors for 4point functions" << std::endl;
        exit(-1);
      }
    }
  }

  // Build all combinations of random vectors for 4point functions.
  // if two or more quarks have the same id, skip combinations where at least
  // two are equal
  for(size_t rnd1 = 0; rnd1 < rndq1; ++rnd1) {
    for(size_t rnd2 = 0; rnd2 < rndq2; ++rnd2) {
      if( (q1 != q2 ) || ( (q1 == q2) && (rnd1 != rnd2) ) ) {
        for(size_t rnd3 = 0; rnd3 < rndq3; ++rnd3) {
          if( ( (q1 != q3) && (q2 != q3) ) ||
              ( (q1 == q3) && (q2 == q3) && (rnd1 != rnd3) && (rnd2 != rnd3) ) ||
              ( (q1 == q3) && (q2 != q3) && (rnd1 != rnd3) ) ||
              ( (q1 != q3) && (q2 == q3) && (rnd2 != rnd3) ) ) {
            for(size_t rnd4 = 0; rnd4 < rndq4; ++rnd4) {
              if( ( (q1 != q4) && (q2 != q4) && (q3 != q4) ) ||
                  ( (q1 == q4) && (q2 == q4) && (q3 == q4) && (rnd1 != rnd4) && (rnd2 != rnd4) && (rnd3 != rnd4) ) ||
                  ( (q1 == q4) && (q2 != q4) && (q3 != q4) && (rnd1 != rnd4) ) ||
                  ( (q1 != q4) && (q2 == q4) && (q3 != q4) && (rnd2 != rnd4) ) ||
                  ( (q1 != q4) && (q2 != q4) && (q3 == q4) && (rnd3 != rnd4) ) ||
                  ( (q1 == q4) && (q2 == q4) && (q3 != q4) && (rnd1 != rnd4) && (rnd2 != rnd4) ) ||
                  ( (q1 == q4) && (q2 != q4) && (q3 == q4) && (rnd1 != rnd4) && (rnd3 != rnd4) ) ||
                  ( (q1 != q4) && (q2 == q4) && (q3 == q4) && (rnd2 != rnd4) && (rnd3 != rnd4) ) ) {
                rnd_vec_4pt.emplace_back(std::array<size_t, 4> {{rnd1, rnd2, rnd3, rnd4}});
              }
            }
          }
        }
      }
    }
  }

//  std::cout << "rnd_vec test: " << rnd_vec_C4.size() << std::endl;
//  for(auto& r : rnd_vec_C4) {
//    std::cout << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << std::endl;
//  }
}

} // end of unnamed namespace

void GlobalData::init_lookup_tables() {

  try{
    init_lookup_corr(correlator_list, operator_list, lookup_corr, lookup_vdv, 
                     lookup_rvdvr, index_of_unity);    
    init_lookup_2pt(correlator_list, operator_list, lookup_corr, lookup_2pt);
    init_lookup_4pt(correlator_list, operator_list, lookup_corr, lookup_4pt,
                    lookup_4pt_3_IO);

    init_lookup_C2plus_IO(correlator_list, operator_list, lookup_corr, lookup_2pt,
                    lookup_2pt_IO);
    init_lookup_C4I2plus_IO(correlator_list, operator_list, lookup_2pt, lookup_4pt,
                    lookup_4pt_1_IO, lookup_4pt_2_IO);

    set_rnd_vec_1pt(quarks, rnd_vec_1pt);
    set_rnd_vec_2pt(quarks, rnd_vec_2pt);
    set_rnd_vec_3pt(quarks, rnd_vec_3pt);
    set_rnd_vec_4pt(quarks, rnd_vec_4pt);

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

