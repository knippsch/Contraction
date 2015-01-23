/*
 * operators.h
 *
 *  Created on: Jan, 2015
 *      Author: knippsch
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <typedefs.h>

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
/// @brief quark type that contains all quark propagator informations
struct quark {

	std::string type;
	int number_of_rnd_vec;
	std::string dilution_T;
	int number_of_dilution_T;
	std::string dilution_E;
	int number_of_dilution_E;
	std::string dilution_D;
	int number_of_dilution_D;

	/// @brief Constructor.
	quark (std::string type, int number_of_rnd_vec, std::string dilution_T,
			int number_of_dilution_T, std::string dilution_E,
			int number_of_dilution_E, std::string dilution_D,
			int number_of_dilution_D) :
			    type(type), number_of_rnd_vec(number_of_rnd_vec), 
          dilution_T(dilution_T), number_of_dilution_T(number_of_dilution_T), 
          dilution_E(dilution_E), number_of_dilution_E(number_of_dilution_E), 
          dilution_D(dilution_D), number_of_dilution_D(number_of_dilution_D) {}
};
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
/// @brief operator type that contains all operator informations
struct Operators {

public: // TODO: should be changed to private at a later point
  std::vector<int> gammas;
  std::array<int, 3> dil_vec;
  std::vector< std::vector<std::array<int, 3> > > mom_vec;

public:
	/// @brief Constructor.
	Operators (std::vector<int> gammas, std::array<int, 3> dil_vec, 
             std::vector<std::vector<std::array<int, 3> > > mom_vec) :
                           gammas(gammas), dil_vec(dil_vec), mom_vec(mom_vec) {}

};
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
/// @brief correlator type that contains all correlator informations
struct Correlators {

public: // TODO: should be changed to private at a later point
  std::string type;
  std::vector<int> quark_numbers;
  std::vector<int> operator_numbers;
  std::string GEVP;
  std::vector<int> tot_mom;

public:
	/// @brief Constructor.
	Correlators(std::string type, std::vector<int> quark_numbers, 
              std::vector<int> operator_numbers, std::string GEVP, 
              std::vector<int> tot_mom) :
                    type(type), quark_numbers(quark_numbers), 
                    operator_numbers(operator_numbers), GEVP(GEVP), 
                    tot_mom(tot_mom) {}

};
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
typedef std::vector<Operators> Operator_list;
typedef std::vector<Correlators> Correlator_list;

#endif /* OPERATORS_H_ */
