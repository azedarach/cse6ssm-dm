// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef SFERMIONS_H
#define SFERMIONS_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace sfermions {

enum Sparticles {
   up = 0,
   down = 1,
   neutrino = 2,
   electron = 3,
   NUMBER_OF_MSSM_SPARTICLES
};

extern const double Isospin[NUMBER_OF_MSSM_SPARTICLES];
extern const double Hypercharge_left[NUMBER_OF_MSSM_SPARTICLES];
extern const double Hypercharge_right[NUMBER_OF_MSSM_SPARTICLES];
extern const double U1prime_charge_left[NUMBER_OF_MSSM_SPARTICLES];
extern const double U1prime_charge_right[NUMBER_OF_MSSM_SPARTICLES];

/**
 * data needed to fill 2 x 2 sfermion mass matrix 
 */ 
struct Mass_data {
   double ml2;    ///< soft mass of left-handed sfermion
   double mr2;    ///< soft mass of right-handed sfermion
   double yf;     ///< Yukawa coupling
   double vd, vu; ///< Higgs VEVs
   double gY, g2; ///< gauge couplings (not GUT normalized)
   double Tyf;    ///< trilinear coupling
   double mu;     ///< Superpotential parameter
   double T3;     ///< weak isospin
   double Yl;     ///< Hypercharge of left-handed sfermion
   double Yr;     ///< Hypercharge of right-handed sfermion
};

/**
 * data needed to fill 2 x 2 sfermion mass matrix, in
 * the case of U(1)-extended models with a single
 * SM singlet charged under the U(1)' (i.e. the E6SSM)
 */
struct U1_extended_mass_data_one_singlet {
   double ml2;        ///< soft mass of left-handed sfermion
   double mr2;        ///< soft mass of right-handed sfermion
   double yf;         ///< Yukawa coupling
   double vd, vu;     ///< Higgs VEVs
   double vs;         ///< SM singlet VEV
   double QHd, QHu;   ///> U(1)'-charges of Higgs doublets
   double QS;         ///< U(1)'-charge of SM singlet
   double gY, g2, gN; ///< gauge couplings (not GUT normalized)
   double Tyf;        ///< trilinear coupling
   double mueff;      ///< superpotential parameter
   double T3;         ///< weak isospin
   double Yl;         ///< Hypercharge of left-handed sfermion
   double Yr;         ///< Hypercharge of right-handed sfermion
   double Ql;         ///< U(1)'-charge of left-handed sfermion
   double Qr;         ///< U(1)'-charge of right-handed sfermion
};

/**
 * data needed to fill 2 x 2 sfermion mass matrix, in
 * the case of U(1)-extended models with a two
 * SM singlets charged under the U(1)' (i.e. the NE6SSM)
 */
struct U1_extended_mass_data_two_singlets {
   double ml2;         ///< soft mass of left-handed sfermion
   double mr2;         ///< soft mass of right-handed sfermion
   double yf;          ///< Yukawa coupling
   double vd, vu;      ///< Higgs VEVs
   double vs, vsb;     ///< SM singlet VEV
   double QHd, QHu;   ///> U(1)'-charges of Higgs doublets
   double QS, QSb;     ///< U(1)'-charges of SM singlets
   double gY, g2, gN;  ///< gauge couplings (not GUT normalized)
   double Tyf;         ///< trilinear coupling
   double mueff;       ///< superpotential parameter
   double T3;          ///< weak isospin
   double Yl;          ///< Hypercharge of left-handed sfermion
   double Yr;          ///< Hypercharge of right-handed sfermion
   double Ql;          ///< U(1)'-charge of left-handed sfermion
   double Qr;          ///< U(1)'-charge of right-handed sfermion
};

double diagonalize_sfermions_2x2(const Mass_data&,
                                 Eigen::Array<double,2,1>&);

double diagonalize_sfermions_2x2(const U1_extended_mass_data_one_singlet&,
                                 Eigen::Array<double,2,1>&);

double diagonalize_sfermions_2x2(const U1_extended_mass_data_two_singlets&,
                                 Eigen::Array<double,2,1>&);

} // namespace sfermions
} // namespace flexiblesusy

#endif
