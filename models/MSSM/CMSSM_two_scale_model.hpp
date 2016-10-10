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

// File generated at Fri 10 Jul 2015 12:05:10

/**
 * @file CMSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Fri 10 Jul 2015 12:05:10 with FlexibleSUSY
 * 1.1.0 (git commit: v1.1.0) and SARAH 4.5.8 .
 */

#ifndef CMSSM_TWO_SCALE_H
#define CMSSM_TWO_SCALE_H

#include "CMSSM_model.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "CMSSM_two_scale_input_parameters.hpp"
#include "two_scale_model.hpp"

#include <iosfwd>
#include <string>

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Two_scale;
/**
 * @class CMSSM<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class CMSSM<Two_scale> : public Two_scale_model, public MSSM_mass_eigenstates {
public:
   explicit CMSSM(const CMSSM_input_parameters<Two_scale>& input_ = CMSSM_input_parameters<Two_scale>());
   virtual ~CMSSM();

   virtual void clear();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(unsigned);
   void set_number_of_ewsb_iterations(std::size_t);
   const CMSSM_input_parameters<Two_scale>& get_input() const;
   void set_input_parameters(const CMSSM_input_parameters<Two_scale>&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level
   double get_ewsb_output_parameter(unsigned) const;

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

   // model interface
   double get_parameter(unsigned) const;
   void set_parameter(unsigned, double);

private:
   struct EWSB_args {
      CMSSM<Two_scale>* model;
      unsigned ewsb_loop_order;
   };

   // input parameters
   CMSSM_input_parameters<Two_scale> input;

   std::size_t number_of_ewsb_iterations;
   unsigned ewsb_loop_order;
   double precision;              ///< RG running precision
   double ewsb_iteration_precision;
   static const std::size_t number_of_ewsb_equations = 2;

   Eigen::Array<double,number_of_tadpole_equations,1> ewsb_solution;

   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(unsigned);
   int solve_ewsb_iteratively_with(EWSB_solver*, const double[number_of_ewsb_equations]);
   int check_ewsb_solution(double);
   void ewsb_initial_guess(double[number_of_ewsb_equations]);
   int ewsb_step(double[number_of_ewsb_equations]) const;
   static int ewsb_step(const gsl_vector*, void*, gsl_vector*);
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   static int tadpole_equations(const gsl_vector*, void*, gsl_vector*);

};

std::ostream& operator<<(std::ostream&, const CMSSM<Two_scale>&);

} // namespace flexiblesusy

#endif
