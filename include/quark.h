/*
 * quark.h
 *
 *  Created on: Jun 20, 2013
 *      Author: knippsch
 */

#ifndef QUARK_H_
#define QUARK_H_

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
			type(type), number_of_rnd_vec(number_of_rnd_vec), dilution_T(dilution_T), number_of_dilution_T(
					number_of_dilution_T), dilution_E(dilution_E), number_of_dilution_E(
					number_of_dilution_E), dilution_D(dilution_D), number_of_dilution_D(
					number_of_dilution_D) {
	}
	/// @brief Destructor.
	virtual ~quark () {
	}
};

#endif /* QUARK_H_ */
