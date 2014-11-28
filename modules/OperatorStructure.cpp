#include "OperatorStructure.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

void LapH::init_from_infile(std::vector<pdg>& op, std::vector<pdg_C2>& op_C2) {

  //TODO: think about the numbers of momenta, displacements and gamma
  const size_t nb_mom = global_data->get_number_of_momenta();
  const size_t nb_mom_sq = global_data->get_number_of_max_mom() + 1;
  //TODO: include displacement into dg (displacementgamma) multiindex
  std::vector<size_t> dg {0, 5};
  const size_t nb_dg = dg.size();

  // nb_op - number of combinations of three-momenta and gamma structures
  // op    - vector of all three-momenta, three-displacements and gamma 
  //         structure combinations
  const size_t nb_op = nb_mom * nb_dg;
  op.resize(nb_op);
  set_default(op);

  // nb_op_C2 - number of combinations of absolute values squared of momenta
  //            and gamma-displacement combinations for 2pt-fct
  // op_C2    - vector of all combinations for 2pt-fct and vector of 
  //            op-index-pairs with all corresponding three-vectors and gammas
  const size_t nb_op_C2 = nb_mom_sq * nb_dg * nb_dg;
  op_C2.resize(nb_op_C2);
  set_default(op, op_C2);

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

// function to initialize the vector of all necessary operators
void LapH::set_default(std::vector<pdg>& op){

  const int max_mom_squared = global_data->get_number_of_max_mom();
  const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();
  std::vector<size_t> dg {0, 5};
  const size_t nb_dg = dg.size();

  size_t i = 0;

  // all three-momenta
  for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
    for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
      for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
        if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
          continue;
        }
        // all gamma structures
        for(size_t gam = 0; gam < nb_dg; gam++){

          op[i].p[0] = ipx;
          op[i].p[1] = ipy;
          op[i].p[2] = ipz;
          op[i].dis[0] = 0;
          op[i].dis[1] = 0;
          op[i].dis[2] = 0;
          op[i].gamma[0] = dg[gam];
          op[i].gamma[1] = 4;
          op[i].gamma[2] = 4;
          op[i].gamma[3] = 4;
  
          op[i].id = i;
          i++;
        }
      }
    }
  }

  if(i != op.size()){
    std::cout << "Error in LapH::set_default(): nb_op not equal to allocated "
                 "number of operators" << std::endl;
    exit(0);
  } 

}

// function to obtain the index combinations in op for the 2pt-fct
void LapH::set_default(std::vector<pdg>& op, std::vector<pdg_C2>& op_C2){

  const size_t nb_mom_sq = global_data->get_number_of_max_mom() + 1;
  std::vector<size_t> dg {0, 5};
  const size_t nb_dg = dg.size();

  size_t nb_op = op.size();
  size_t j = 0;

  // run over all momenta squared (back-to-back hardcoded) and gamma 
  // combinations
  for(size_t p_sq = 0; p_sq < nb_mom_sq; p_sq++){
    for(size_t so = 0; so < nb_dg; so++){
      for(size_t si = 0; si < nb_dg; si++){

        // index for access of element
        size_t i = p_sq * nb_dg*nb_dg + so * nb_dg + si;

        // save p^2 and gamma structure at source and sink
        op_C2[i].p_sq = p_sq;
        op_C2[i].dg_so = dg[so];
        op_C2[i].dg_si = dg[si];

        // loop over op and set index pairs
        for(auto& el : op)
          if((el.p[0]*el.p[0] + el.p[1]*el.p[1] + el.p[2]*el.p[2]) == p_sq){ 
            if(el.gamma[0] == dg[so]){
              size_t id1 = el.id;
              // thats the generalized version of nb_mom - p - 1 including 
              // a faster running gamma structure
              size_t id2 = nb_op - nb_dg * (id1/nb_dg + 1) + si;
              // warning because array has no list-of-arrays constructor but 
              // works. Can change this to pair structure.
              op_C2[i].index.emplace_back(std::pair<size_t, size_t>(id1, id2));
            }
          }

        j++;
        
      }
    }
  }

  if(j != op_C2.size()){
    std::cout << "Error in LapH::set_default(): nb_op not equal to allocated "
                 "number of operators" << std::endl;
    exit(0);
  } 

}




