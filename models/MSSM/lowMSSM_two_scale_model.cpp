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

// File generated at Wed 15 Jul 2015 13:12:37

/**
 * @file lowMSSM_two_scale_model.cpp
 * @brief implementation of the lowMSSM model class
 *
 * Contains the definition of the lowMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 15 Jul 2015 13:12:37 with FlexibleSUSY
 * 1.1.0 (git commit: v1.1.0) and SARAH 4.5.8 .
 */

#include "lowMSSM_two_scale_model.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "functors.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace MSSM_info;

#define CLASSNAME lowMSSM<Two_scale>

CLASSNAME::lowMSSM(const lowMSSM_input_parameters<Two_scale>& input_)
   : Two_scale_model()
   , MSSM_mass_eigenstates()
   , input(input_)
   , number_of_ewsb_iterations(200)
   , ewsb_loop_order(2)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , ewsb_solution(Eigen::Array<double,number_of_tadpole_equations,1>::Zero())

{
}

CLASSNAME::~lowMSSM()
{
}

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
   diagonalization_precision = precision_;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

double CLASSNAME::get_ewsb_output_parameter(unsigned i) const
{
   return ewsb_solution(i);
}

const lowMSSM_input_parameters<Two_scale>& CLASSNAME::get_input() const
{
   return input;
}

void CLASSNAME::set_input_parameters(const lowMSSM_input_parameters<Two_scale>& input_)
{
   input = input_;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[2];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] -= two_loop_tadpole[0];
         tadpole[1] -= two_loop_tadpole[1];

      }
   }
}

/**
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::EWSB_args struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::EWSB_args struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   lowMSSM<Two_scale>* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mHd2(gsl_vector_get(x, 0));
   model->set_mHu2(gsl_vector_get(x, 1));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double tadpole[number_of_ewsb_equations] = { 0. };

   model->tadpole_equations(tadpole);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return is_finite<number_of_ewsb_equations>(tadpole) ? GSL_SUCCESS : GSL_EDOM;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   EWSB_solver* solvers[] = {
      new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrid),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden)
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   double x_init[number_of_ewsb_equations] = { 0. };
   ewsb_initial_guess(x_init);

#ifdef ENABLE_VERBOSE
   std::cout << "Solving EWSB equations ...\n"
      "\tInitial guess: x_init =";
   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      std::cout << ' ' << x_init[i];
   std::cout << '\n';
#endif

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i], x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   for_each(solvers, solvers + number_of_solvers, Delete_object());

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const double x_init[number_of_ewsb_equations]
)
{
   const int status = solver->solve(x_init);

   mHd2 = solver->get_solution(0);
   mHu2 = solver->get_solution(1);


   return status;
}

int CLASSNAME::check_ewsb_solution(double precision)
{
   double tadpole[number_of_ewsb_equations];

   if (ewsb_loop_order > 0)
      calculate_DRbar_masses();

   tadpole_equations(tadpole);

   double residual = Abs(tadpole[0]) + Abs(tadpole[1]);

   return (residual < precision ? EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) - 3*Power(
      vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*
      Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) - 3*Power(
      vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*
      Sqr(vd)))/vu);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

void CLASSNAME::ewsb_initial_guess(double x_init[number_of_ewsb_equations])
{
   x_init[0] = mHd2;
   x_init[1] = mHu2;

}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param ewsb_parameters new EWSB output parameters.  \a
 * ewsb_parameters is only modified if all new parameters are finite.
 *
 * @return GSL_SUCCESS if new EWSB output parameters are finite,
 * GSL_EDOM otherwise.
 */
int CLASSNAME::ewsb_step(double ewsb_parameters[number_of_ewsb_equations]) const
{
   int error;
   double tadpole[number_of_ewsb_equations] = { 0. };

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh(0));
      tadpole[1] += Re(tadpole_hh(1));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[2];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] += two_loop_tadpole[0];
         tadpole[1] += two_loop_tadpole[1];

      }
   }

   double mHd2;
   double mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) + 40*
      tadpole[0] - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) + 3*vd*Sqr(g1)*
      Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) + 40*
      tadpole[1] - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*
      Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd)))/vu);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);


   if (is_finite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = mHd2;
      ewsb_parameters[1] = mHu2;

   } else {
      error = GSL_EDOM;
   }

   return error;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param x old EWSB output parameters
 * @param params further function parameters
 * @param f new EWSB output parameters
 *
 * @return Returns status of CLASSNAME::ewsb_step
 */
int CLASSNAME::ewsb_step(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   lowMSSM<Two_scale>* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mHd2 = gsl_vector_get(x, 0);
   const double mHu2 = gsl_vector_get(x, 1);

   model->set_mHd2(mHd2);
   model->set_mHu2(mHu2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { mHd2, mHu2 };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "lowMSSM (solver type: two_scale)\n"
           "========================================\n";
   MSSM_mass_eigenstates::print(ostr);
}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   MSSM_mass_eigenstates::clear();
}

std::string CLASSNAME::name() const
{
   return "lowMSSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   MSSM_mass_eigenstates::run_to(scale, eps);
}

double CLASSNAME::get_parameter(unsigned parameter) const
{
   if (parameter >= MSSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case MSSM_info::Yd00:
      return Yd(0,0);
   case MSSM_info::Yd01:
      return Yd(0,1);
   case MSSM_info::Yd02:
      return Yd(0,2);
   case MSSM_info::Yd10:
      return Yd(1,0);
   case MSSM_info::Yd11:
      return Yd(1,1);
   case MSSM_info::Yd12:
      return Yd(1,2);
   case MSSM_info::Yd20:
      return Yd(2,0);
   case MSSM_info::Yd21:
      return Yd(2,1);
   case MSSM_info::Yd22:
      return Yd(2,2);
   case MSSM_info::Ye00:
      return Ye(0,0);
   case MSSM_info::Ye01:
      return Ye(0,1);
   case MSSM_info::Ye02:
      return Ye(0,2);
   case MSSM_info::Ye10:
      return Ye(1,0);
   case MSSM_info::Ye11:
      return Ye(1,1);
   case MSSM_info::Ye12:
      return Ye(1,2);
   case MSSM_info::Ye20:
      return Ye(2,0);
   case MSSM_info::Ye21:
      return Ye(2,1);
   case MSSM_info::Ye22:
      return Ye(2,2);
   case MSSM_info::Yu00:
      return Yu(0,0);
   case MSSM_info::Yu01:
      return Yu(0,1);
   case MSSM_info::Yu02:
      return Yu(0,2);
   case MSSM_info::Yu10:
      return Yu(1,0);
   case MSSM_info::Yu11:
      return Yu(1,1);
   case MSSM_info::Yu12:
      return Yu(1,2);
   case MSSM_info::Yu20:
      return Yu(2,0);
   case MSSM_info::Yu21:
      return Yu(2,1);
   case MSSM_info::Yu22:
      return Yu(2,2);
   case MSSM_info::Mu:
      return Mu;
   case MSSM_info::g1:
      return g1;
   case MSSM_info::g2:
      return g2;
   case MSSM_info::g3:
      return g3;
   case MSSM_info::vd:
      return vd;
   case MSSM_info::vu:
      return vu;
   case MSSM_info::TYd00:
      return TYd(0,0);
   case MSSM_info::TYd01:
      return TYd(0,1);
   case MSSM_info::TYd02:
      return TYd(0,2);
   case MSSM_info::TYd10:
      return TYd(1,0);
   case MSSM_info::TYd11:
      return TYd(1,1);
   case MSSM_info::TYd12:
      return TYd(1,2);
   case MSSM_info::TYd20:
      return TYd(2,0);
   case MSSM_info::TYd21:
      return TYd(2,1);
   case MSSM_info::TYd22:
      return TYd(2,2);
   case MSSM_info::TYe00:
      return TYe(0,0);
   case MSSM_info::TYe01:
      return TYe(0,1);
   case MSSM_info::TYe02:
      return TYe(0,2);
   case MSSM_info::TYe10:
      return TYe(1,0);
   case MSSM_info::TYe11:
      return TYe(1,1);
   case MSSM_info::TYe12:
      return TYe(1,2);
   case MSSM_info::TYe20:
      return TYe(2,0);
   case MSSM_info::TYe21:
      return TYe(2,1);
   case MSSM_info::TYe22:
      return TYe(2,2);
   case MSSM_info::TYu00:
      return TYu(0,0);
   case MSSM_info::TYu01:
      return TYu(0,1);
   case MSSM_info::TYu02:
      return TYu(0,2);
   case MSSM_info::TYu10:
      return TYu(1,0);
   case MSSM_info::TYu11:
      return TYu(1,1);
   case MSSM_info::TYu12:
      return TYu(1,2);
   case MSSM_info::TYu20:
      return TYu(2,0);
   case MSSM_info::TYu21:
      return TYu(2,1);
   case MSSM_info::TYu22:
      return TYu(2,2);
   case MSSM_info::BMu:
      return BMu;
   case MSSM_info::mq200:
      return mq2(0,0);
   case MSSM_info::mq201:
      return mq2(0,1); 
   case MSSM_info::mq202:
      return mq2(0,2); 
   case MSSM_info::mq210:
      return mq2(1,0);
   case MSSM_info::mq211:
      return mq2(1,1); 
   case MSSM_info::mq212:
      return mq2(1,2);
   case MSSM_info::mq220:
      return mq2(2,0);
   case MSSM_info::mq221:
      return mq2(2,1);
   case MSSM_info::mq222:
      return mq2(2,2); 
   case MSSM_info::ml200:
      return ml2(0,0);
   case MSSM_info::ml201:
      return ml2(0,1); 
   case MSSM_info::ml202:
      return ml2(0,2); 
   case MSSM_info::ml210:
      return ml2(1,0);
   case MSSM_info::ml211:
      return ml2(1,1); 
   case MSSM_info::ml212:
      return ml2(1,2);
   case MSSM_info::ml220:
      return ml2(2,0);
   case MSSM_info::ml221:
      return ml2(2,1);
   case MSSM_info::ml222:
      return ml2(2,2); 
   case MSSM_info::mHd2:
      return mHd2;
   case MSSM_info::mHu2:
      return mHu2;
   case MSSM_info::md200:
      return md2(0,0);
   case MSSM_info::md201:
      return md2(0,1); 
   case MSSM_info::md202:
      return md2(0,2); 
   case MSSM_info::md210:
      return md2(1,0);
   case MSSM_info::md211:
      return md2(1,1); 
   case MSSM_info::md212:
      return md2(1,2);
   case MSSM_info::md220:
      return md2(2,0);
   case MSSM_info::md221:
      return md2(2,1);
   case MSSM_info::md222:
      return md2(2,2); 
   case MSSM_info::mu200:
      return mu2(0,0);
   case MSSM_info::mu201:
      return mu2(0,1); 
   case MSSM_info::mu202:
      return mu2(0,2); 
   case MSSM_info::mu210:
      return mu2(1,0);
   case MSSM_info::mu211:
      return mu2(1,1); 
   case MSSM_info::mu212:
      return mu2(1,2);
   case MSSM_info::mu220:
      return mu2(2,0);
   case MSSM_info::mu221:
      return mu2(2,1);
   case MSSM_info::mu222:
      return mu2(2,2); 
   case MSSM_info::me200:
      return me2(0,0);
   case MSSM_info::me201:
      return me2(0,1); 
   case MSSM_info::me202:
      return me2(0,2); 
   case MSSM_info::me210:
      return me2(1,0);
   case MSSM_info::me211:
      return me2(1,1); 
   case MSSM_info::me212:
      return me2(1,2);
   case MSSM_info::me220:
      return me2(2,0);
   case MSSM_info::me221:
      return me2(2,1);
   case MSSM_info::me222:
      return me2(2,2); 
   case MSSM_info::MassB:
      return MassB;
   case MSSM_info::MassWB:
      return MassWB;
   case MSSM_info::MassG:
      return MassG;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

void CLASSNAME::set_parameter(unsigned parameter, double x)
{
   if (parameter >= MSSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case MSSM_info::Yd00:
      Yd(0,0) = x;
      break;
   case MSSM_info::Yd01:
      Yd(0,1) = x;
      break;
   case MSSM_info::Yd02:
      Yd(0,2) = x;
      break;
   case MSSM_info::Yd10:
      Yd(1,0) = x;
      break;
   case MSSM_info::Yd11:
      Yd(1,1) = x;
      break;
   case MSSM_info::Yd12:
      Yd(1,2) = x;
      break;
   case MSSM_info::Yd20:
      Yd(2,0) = x;
      break;
   case MSSM_info::Yd21:
      Yd(2,1) = x;
      break;
   case MSSM_info::Yd22:
      Yd(2,2) = x;
      break;
   case MSSM_info::Ye00:
      Ye(0,0) = x;
      break;
   case MSSM_info::Ye01:
      Ye(0,1) = x;
      break;
   case MSSM_info::Ye02:
      Ye(0,2) = x;
      break;
   case MSSM_info::Ye10:
      Ye(1,0) = x;
      break;
   case MSSM_info::Ye11:
      Ye(1,1) = x;
      break;
   case MSSM_info::Ye12:
      Ye(1,2) = x;
      break;
   case MSSM_info::Ye20:
      Ye(2,0) = x;
      break;
   case MSSM_info::Ye21:
      Ye(2,1) = x;
      break;
   case MSSM_info::Ye22:
      Ye(2,2) = x;
      break;
   case MSSM_info::Yu00:
      Yu(0,0) = x;
      break;
   case MSSM_info::Yu01:
      Yu(0,1) = x;
      break;
   case MSSM_info::Yu02:
      Yu(0,2) = x;
      break;
   case MSSM_info::Yu10:
      Yu(1,0) = x;
      break;
   case MSSM_info::Yu11:
      Yu(1,1) = x;
      break;
   case MSSM_info::Yu12:
      Yu(1,2) = x;
      break;
   case MSSM_info::Yu20:
      Yu(2,0) = x;
      break;
   case MSSM_info::Yu21:
      Yu(2,1) = x;
      break;
   case MSSM_info::Yu22:
      Yu(2,2) = x;
      break;
   case MSSM_info::Mu:
      Mu = x;
      break;
   case MSSM_info::g1:
      g1 = x;
      break;
   case MSSM_info::g2:
      g2 = x;
      break;
   case MSSM_info::g3:
      g3 = x;
      break;
   case MSSM_info::vd:
      vd = x;
      break;
   case MSSM_info::vu:
      vu = x;
      break;
   case MSSM_info::TYd00:
      TYd(0,0) = x;
      break;
   case MSSM_info::TYd01:
      TYd(0,1) = x;
      break;
   case MSSM_info::TYd02:
      TYd(0,2) = x;
      break;
   case MSSM_info::TYd10:
      TYd(1,0) = x;
      break;
   case MSSM_info::TYd11:
      TYd(1,1) = x;
      break;
   case MSSM_info::TYd12:
      TYd(1,2) = x;
      break;
   case MSSM_info::TYd20:
      TYd(2,0) = x;
      break;
   case MSSM_info::TYd21:
      TYd(2,1) = x;
      break;
   case MSSM_info::TYd22:
      TYd(2,2) = x;
      break;
   case MSSM_info::TYe00:
      TYe(0,0) = x;
      break;
   case MSSM_info::TYe01:
      TYe(0,1) = x;
      break;
   case MSSM_info::TYe02:
      TYe(0,2) = x;
      break;
   case MSSM_info::TYe10:
      TYe(1,0) = x;
      break;
   case MSSM_info::TYe11:
      TYe(1,1) = x;
      break;
   case MSSM_info::TYe12:
      TYe(1,2) = x;
      break;
   case MSSM_info::TYe20:
      TYe(2,0) = x;
      break;
   case MSSM_info::TYe21:
      TYe(2,1) = x;
      break;
   case MSSM_info::TYe22:
      TYe(2,2) = x;
      break;
   case MSSM_info::TYu00:
      TYu(0,0) = x;
      break;
   case MSSM_info::TYu01:
      TYu(0,1) = x;
      break;
   case MSSM_info::TYu02:
      TYu(0,2) = x;
      break;
   case MSSM_info::TYu10:
      TYu(1,0) = x;
      break;
   case MSSM_info::TYu11:
      TYu(1,1) = x;
      break;
   case MSSM_info::TYu12:
      TYu(1,2) = x;
      break;
   case MSSM_info::TYu20:
      TYu(2,0) = x;
      break;
   case MSSM_info::TYu21:
      TYu(2,1) = x;
      break;
   case MSSM_info::TYu22:
      TYu(2,2) = x;
      break;
   case MSSM_info::BMu:
      BMu = x;
      break;
   case MSSM_info::mq200:
      mq2(0,0) = x;
      break;
   case MSSM_info::mq201:
      mq2(0,1) = x;
      break; 
   case MSSM_info::mq202:
      mq2(0,2) = x;
      break; 
   case MSSM_info::mq210:
      mq2(1,0) = x;
      break;
   case MSSM_info::mq211:
      mq2(1,1) = x;
      break; 
   case MSSM_info::mq212:
      mq2(1,2) = x;
      break;
   case MSSM_info::mq220:
      mq2(2,0) = x;
      break;
   case MSSM_info::mq221:
      mq2(2,1) = x;
      break;
   case MSSM_info::mq222:
      mq2(2,2) = x;
      break; 
   case MSSM_info::ml200:
      ml2(0,0) = x;
      break;
   case MSSM_info::ml201:
      ml2(0,1) = x;
      break; 
   case MSSM_info::ml202:
      ml2(0,2) = x;
      break; 
   case MSSM_info::ml210:
      ml2(1,0) = x;
      break;
   case MSSM_info::ml211:
      ml2(1,1) = x;
      break;
   case MSSM_info::ml212:
      ml2(1,2) = x;
      break;
   case MSSM_info::ml220:
      ml2(2,0) = x;
      break;
   case MSSM_info::ml221:
      ml2(2,1) = x;
      break;
   case MSSM_info::ml222:
      ml2(2,2) = x;
      break; 
   case MSSM_info::mHd2:
      mHd2 = x;
      break;
   case MSSM_info::mHu2:
      mHu2 = x;
      break;
   case MSSM_info::md200:
      md2(0,0) = x;
      break;
   case MSSM_info::md201:
      md2(0,1) = x;
      break; 
   case MSSM_info::md202:
      md2(0,2) = x;
      break; 
   case MSSM_info::md210:
      md2(1,0) = x;
      break;
   case MSSM_info::md211:
      md2(1,1) = x;
      break; 
   case MSSM_info::md212:
      md2(1,2) = x;
      break;
   case MSSM_info::md220:
      md2(2,0) = x;
      break;
   case MSSM_info::md221:
      md2(2,1) = x;
      break;
   case MSSM_info::md222:
      md2(2,2) = x;
      break; 
   case MSSM_info::mu200:
      mu2(0,0) = x;
      break;
   case MSSM_info::mu201:
      mu2(0,1) = x;
      break; 
   case MSSM_info::mu202:
      mu2(0,2) = x;
      break; 
   case MSSM_info::mu210:
      mu2(1,0) = x;
      break;
   case MSSM_info::mu211:
      mu2(1,1) = x;
      break; 
   case MSSM_info::mu212:
      mu2(1,2) = x;
      break;
   case MSSM_info::mu220:
      mu2(2,0) = x;
      break;
   case MSSM_info::mu221:
      mu2(2,1) = x;
      break;
   case MSSM_info::mu222:
      mu2(2,2) = x;
      break; 
   case MSSM_info::me200:
      me2(0,0) = x;
      break;
   case MSSM_info::me201:
      me2(0,1) = x;
      break; 
   case MSSM_info::me202:
      me2(0,2) = x;
      break; 
   case MSSM_info::me210:
      me2(1,0) = x;
      break;
   case MSSM_info::me211:
      me2(1,1) = x;
      break; 
   case MSSM_info::me212:
      me2(1,2) = x;
      break;
   case MSSM_info::me220:
      me2(2,0) = x;
      break;
   case MSSM_info::me221:
      me2(2,1) = x;
      break;
   case MSSM_info::me222:
      me2(2,2) = x;
      break; 
   case MSSM_info::MassB:
      MassB = x;
      break;
   case MSSM_info::MassWB:
      MassWB = x;
      break;
   case MSSM_info::MassG:
      MassG = x;
      break;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

std::ostream& operator<<(std::ostream& ostr, const lowMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
