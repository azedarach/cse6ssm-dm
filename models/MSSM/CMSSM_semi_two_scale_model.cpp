
/**
 * @file CMSSM_semi_two_scale_model.cpp
 * @brief implementation of the CMSSM semianalytic model class
 *
 * Contains the definition of the CMSSM semianalytic model class
 * methods which solve EWSB and calculate pole masses and mixings
 * from DRbar parameters.
 */

#include "CMSSM_semi_two_scale_model.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "derivs_root_finder.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "functors.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

#include <iomanip>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace MSSM_info;

#define CLASSNAME CMSSM_semianalytic<Two_scale>

#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

CLASSNAME::CMSSM_semianalytic(const CMSSM_semianalytic_input_parameters<Two_scale>& input_)
   : Two_scale_model()
   , MSSM_mass_eigenstates()
   , input(input_)
   , number_of_ewsb_iterations(200)
   , ewsb_loop_order(2)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , ewsb_solution(Eigen::Array<double,number_of_tadpole_equations,1>::Zero())
   , TYd_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYd_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYe_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYe_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYu_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYu_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , BMu_BMu0_coeff(0), BMu_Azero_coeff(0), BMu_m12_coeff(0)
   , mq2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mHd2_m02_coeff(0), mHd2_m122_coeff(0), mHd2_Azerom12_coeff(0)
   , mHd2_Azero2_coeff(0), mHu2_m02_coeff(0), mHu2_m122_coeff(0)
   , mHu2_Azerom12_coeff(0), mHu2_Azero2_coeff(0)
   , md2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , MassB_Azero_coeff(0), MassB_m12_coeff(0), MassWB_Azero_coeff(0)
   , MassWB_m12_coeff(0), MassG_Azero_coeff(0), MassG_m12_coeff(0)
{
}

CLASSNAME::~CMSSM_semianalytic()
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
   // note: for large masses temporarily use fixed EWSB precision
//   ewsb_iteration_precision = precision_;
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

const CMSSM_semianalytic_input_parameters<Two_scale>& CLASSNAME::get_input() const
{
   return input;
}

CMSSM_semianalytic_input_parameters<Two_scale>& CLASSNAME::get_input()
{
   return input;
}

void CLASSNAME::set_input_parameters(const CMSSM_semianalytic_input_parameters<Two_scale>& input_)
{
   input = input_;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_tadpole_equations]) const
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
   CMSSM_semianalytic<Two_scale>* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   // EWSB output parameters are m0, BMu0
   const double m12 = model->get_input().m12;
   const double Azero = model->get_input().Azero;

   const double m0Sqr = gsl_vector_get(x, 0);
   const double BMu0 = gsl_vector_get(x, 1);

   model->set_soft_parameters_at_current_scale(m0Sqr, m12, Azero, BMu0);

   if (ewsb_loop_order > 0) {
      Problems<MSSM_info::NUMBER_OF_PARTICLES> old_problems = model->get_problems();
      model->calculate_DRbar_masses();
      Problems<MSSM_info::NUMBER_OF_PARTICLES> new_problems = model->get_problems();
      if (new_problems.have_tachyon() && !old_problems.have_tachyon()) {
         for (unsigned i = 0; i < MSSM_info::NUMBER_OF_PARTICLES; ++i) {
            // unflag tachyons introduced during EWSB iteration
            if (new_problems.is_tachyon(i) && !old_problems.is_tachyon(i))
               model->get_problems().unflag_tachyon(i);
         }
      }
   }

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
      has_previous_ewsb_solution = true;
   } else {
      problems.flag_no_ewsb();
      has_previous_ewsb_solution = false;
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
   int error = solver->solve(x_init);

   const double m12 = LOCALINPUT(m12);
   const double Azero = LOCALINPUT(Azero);

   ewsb_solution(0) = solver->get_solution(0);
   ewsb_solution(1) = solver->get_solution(1);

   set_soft_parameters_at_current_scale(ewsb_solution(0), m12, Azero,
      ewsb_solution(1));

   return error;
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

   const double old_m0Sqr = ewsb_solution(0);
   const double old_BMu0 = ewsb_solution(1);

   const double m12 = input.m12;
   const double Azero = input.Azero;

   const double m0Sqr = (AbsSqr(Mu) * (Sqr(vu) - Sqr(vd)) + 0.125 *
      (Sqr(g2) + 0.6 * Sqr(g1)) * (Sqr(vd) + Sqr(vu)) * (Sqr(vu) - Sqr(vd))
      + Sqr(m12) * (mHu2_m122_coeff * Sqr(vu) - mHd2_m122_coeff * Sqr(vd))
      + m12 * Azero * (mHu2_Azerom12_coeff * Sqr(vu) - mHd2_Azerom12_coeff
      * Sqr(vd)) + Sqr(Azero) * (mHu2_Azero2_coeff * Sqr(vu)
      - mHd2_Azero2_coeff * Sqr(vd))) / (mHd2_m02_coeff * Sqr(vd)
      - mHu2_m02_coeff * Sqr(vu));
   const double BMu0 = (-BMu_Azero_coeff * Azero - BMu_m12_coeff * m12
      + vd * vu * (2.0 * AbsSqr(Mu) + mHd2_m02_coeff * m0Sqr +
      mHd2_m122_coeff * Sqr(m12) + mHd2_Azerom12_coeff * m12 * Azero
      + mHd2_Azero2_coeff * Sqr(Azero) + mHu2_m02_coeff * m0Sqr +
      mHu2_m122_coeff * Sqr(m12) + mHu2_Azerom12_coeff * m12 * Azero
      + mHu2_Azero2_coeff * Sqr(Azero)) / (Sqr(vd) + Sqr(vu))) / BMu_BMu0_coeff;

   const bool is_finite = IsFinite(BMu0) && IsFinite(m0Sqr);

   if (!is_finite) {
      ewsb_solution(1) = old_BMu0;
      ewsb_solution(0) = old_m0Sqr;
      error = 1;
   } else {
      ewsb_solution(0) = m0Sqr;
      ewsb_solution(1) = BMu0;
      has_previous_ewsb_solution = true;
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

int CLASSNAME::ewsb_initial_guess(double x_init[number_of_ewsb_equations])
{
   int status = EWSB_solver::SUCCESS;

   x_init[0] = ewsb_solution(0);
   x_init[1] = ewsb_solution(1);

   return status;
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
int CLASSNAME::ewsb_step(double ewsb_parameters[number_of_ewsb_equations])
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

   const double m12 = input.m12;
   const double Azero = input.Azero;

   double m0Sqr;
   double BMu0;

   m0Sqr = (AbsSqr(Mu) * (Sqr(vu) - Sqr(vd)) + 0.125 *
      (Sqr(g2) + 0.6 * Sqr(g1)) * (Sqr(vd) + Sqr(vu)) * (Sqr(vu) - Sqr(vd))
      + Sqr(m12) * (mHu2_m122_coeff * Sqr(vu) - mHd2_m122_coeff * Sqr(vd))
      + m12 * Azero * (mHu2_Azerom12_coeff * Sqr(vu) - mHd2_Azerom12_coeff
      * Sqr(vd)) + Sqr(Azero) * (mHu2_Azero2_coeff * Sqr(vu)
      - mHd2_Azero2_coeff * Sqr(vd)) + tadpole[0] * vd - tadpole[1] * vu) /
      (mHd2_m02_coeff * Sqr(vd) - mHu2_m02_coeff * Sqr(vu));
   BMu0 = (-BMu_Azero_coeff * Azero - BMu_m12_coeff * m12
      + vd * vu * (2.0 * AbsSqr(Mu) + mHd2_m02_coeff * m0Sqr +
      mHd2_m122_coeff * Sqr(m12) + mHd2_Azerom12_coeff * m12 * Azero
      + mHd2_Azero2_coeff * Sqr(Azero) + mHu2_m02_coeff * m0Sqr +
      mHu2_m122_coeff * Sqr(m12) + mHu2_Azerom12_coeff * m12 * Azero
      + mHu2_Azero2_coeff * Sqr(Azero) - tadpole[0] / vd
      - tadpole[1] / vu) / (Sqr(vd) + Sqr(vu))) / BMu_BMu0_coeff;

   const bool is_finite = IsFinite(m0Sqr) && IsFinite(BMu0);

   if (is_finite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = m0Sqr;
      ewsb_parameters[1] = BMu0;
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
   CMSSM_semianalytic* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double m12 = INPUT(m12);
   const double Azero = INPUT(Azero);

   const double m0Sqr = gsl_vector_get(x, 0);
   const double BMu0 = gsl_vector_get(x, 1);

   model->set_soft_parameters_at_current_scale(m0Sqr, m12, Azero, BMu0);

   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { m0Sqr, BMu0 };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "CMSSM_semianalytic\n"
           "(solver type: two_scale)\n"
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
   return "CMSSM_semianalytic";
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

void CLASSNAME::calculate_coefficients(double input_scale)
{
   static const double fit_Azero_values[number_of_fit_points] = {0., 1., 0., 1.};
   static const double fit_m12_values[number_of_fit_points] = {1., 0., 0., 1.};
   static const double fit_m0Sqr_values[number_of_fit_points] = {0., 0., 1., 0.};
   static const double fit_BMu0_values[number_of_fit_points] = {0., 0., 1., 0.};

   // save current set of parameters
   MSSM_soft_parameters saved_pars;
   saved_pars.set(get());

   const double current_scale = get_scale();

   run_to(input_scale, precision);

   MSSM_soft_parameters input_scale_pars;
   input_scale_pars.set(get());

   std::vector<MSSM_soft_parameters> parameter_values;

   Eigen::Matrix<double,number_of_fit_points,2> dimension_one_inputs;
   Eigen::Matrix<double,number_of_fit_points,4> dimension_two_inputs;
   Eigen::Matrix<double,number_of_fit_points,3> soft_bilinear_inputs;

   for (std::size_t i = 0; i < number_of_fit_points; ++i) {
      dimension_one_inputs(i,0) = fit_Azero_values[i];
      dimension_one_inputs(i,1) = fit_m12_values[i];

      dimension_two_inputs(i,0) = Sqr(fit_m0Sqr_values[i]);
      dimension_two_inputs(i,1) = Sqr(fit_m12_values[i]);
      dimension_two_inputs(i,2) = fit_m12_values[i] * fit_Azero_values[i];
      dimension_two_inputs(i,3) = Sqr(fit_Azero_values[i]);

      soft_bilinear_inputs(i,0) = fit_Azero_values[i];
      soft_bilinear_inputs(i,1) = fit_m12_values[i];
      soft_bilinear_inputs(i,2) = fit_BMu0_values[i];

      set_soft_parameters_at_input_scale(fit_m0Sqr_values[i], fit_m12_values[i],
                                         fit_Azero_values[i], fit_BMu0_values[i]);

      run_to(current_scale, precision);

      MSSM_soft_parameters params;
      params.set(get());

      parameter_values.push_back(params);

      set(input_scale_pars.get());
      set_scale(input_scale);
   }

   // solve for coefficients using least squares
   // for implementation in FS, use FS SVD routines if possible
   Eigen::JacobiSVD<Eigen::Matrix<double,number_of_fit_points,2> > dimension_one_svd(dimension_one_inputs, Eigen::ComputeFullU | Eigen::ComputeFullV);
   Eigen::JacobiSVD<Eigen::Matrix<double,number_of_fit_points,4> > dimension_two_svd(dimension_two_inputs, Eigen::ComputeFullU | Eigen::ComputeFullV);
   Eigen::JacobiSVD<Eigen::Matrix<double,number_of_fit_points,3> > soft_bilinear_svd(soft_bilinear_inputs, Eigen::ComputeFullU | Eigen::ComputeFullV);

   Eigen::Matrix<double,number_of_fit_points,1> rhs;
   Eigen::Matrix<double,2,1> dimension_one_solution;
   Eigen::Matrix<double,4,1> dimension_two_solution;
   Eigen::Matrix<double,3,1> soft_bilinear_solution;

   // TODO temporary working (awful!) solution, needs to be replaced...
   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYd(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYd_Azero_coeff(i,j) = dimension_one_solution(0);
         TYd_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYe(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYe_Azero_coeff(i,j) = dimension_one_solution(0);
         TYe_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYu(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYu_Azero_coeff(i,j) = dimension_one_solution(0);
         TYu_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mq2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mq2_m02_coeff(i,j) = dimension_two_solution(0);
         mq2_m122_coeff(i,j) = dimension_two_solution(1);
         mq2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mq2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_ml2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         ml2_m02_coeff(i,j) = dimension_two_solution(0);
         ml2_m122_coeff(i,j) = dimension_two_solution(1);
         ml2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         ml2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHd2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHd2_m02_coeff = dimension_two_solution(0);
   mHd2_m122_coeff = dimension_two_solution(1);
   mHd2_Azerom12_coeff = dimension_two_solution(2);
   mHd2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHu2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHu2_m02_coeff = dimension_two_solution(0);
   mHu2_m122_coeff = dimension_two_solution(1);
   mHu2_Azerom12_coeff = dimension_two_solution(2);
   mHu2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_md2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         md2_m02_coeff(i,j) = dimension_two_solution(0);
         md2_m122_coeff(i,j) = dimension_two_solution(1);
         md2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         md2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mu2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mu2_m02_coeff(i,j) = dimension_two_solution(0);
         mu2_m122_coeff(i,j) = dimension_two_solution(1);
         mu2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mu2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_me2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         me2_m02_coeff(i,j) = dimension_two_solution(0);
         me2_m122_coeff(i,j) = dimension_two_solution(1);
         me2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         me2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassB();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassB_Azero_coeff = dimension_one_solution(0);
   MassB_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassWB();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassWB_Azero_coeff = dimension_one_solution(0);
   MassWB_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassG();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassG_Azero_coeff = dimension_one_solution(0);
   MassG_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_BMu();
   }
   soft_bilinear_solution = soft_bilinear_svd.solve(rhs);
   BMu_Azero_coeff = soft_bilinear_solution(0);
   BMu_m12_coeff = soft_bilinear_solution(1);
   BMu_BMu0_coeff = soft_bilinear_solution(2);

   // reset parameters at initial scale
   set(saved_pars.get());
   set_scale(current_scale);

   set_soft_parameters_at_current_scale(ewsb_solution(0), input.m12,
                                        input.Azero, ewsb_solution(1));
}

void CLASSNAME::set_soft_parameters_at_input_scale(double m0Sqr, double m12, double Azero, double BMu0)
{
   TYd = Azero * Yd;
   TYe = Azero * Ye;
   TYu = Azero * Yu;
   mq2 = m0Sqr * UNITMATRIX(3);
   ml2 = m0Sqr * UNITMATRIX(3);
   mHd2 = m0Sqr;
   mHu2 = m0Sqr;
   md2 = m0Sqr * UNITMATRIX(3);
   mu2 = m0Sqr * UNITMATRIX(3);
   me2 = m0Sqr * UNITMATRIX(3);
   MassB = m12;
   MassWB = m12;
   MassG = m12;
   BMu = BMu0;
}

void CLASSNAME::set_soft_parameters_at_current_scale(double m0Sqr, double m12, double Azero, double BMu0)
{
   TYd = TYd_Azero_coeff * Azero + TYd_m12_coeff * m12;
   TYe = TYe_Azero_coeff * Azero + TYe_m12_coeff * m12;
   TYu = TYu_Azero_coeff * Azero + TYu_m12_coeff * m12;
   BMu = BMu_BMu0_coeff * BMu0 + BMu_Azero_coeff * Azero + BMu_m12_coeff * m12;
   mq2 = mq2_m02_coeff * m0Sqr + mq2_m122_coeff * Sqr(m12)
      + mq2_Azerom12_coeff * Azero * m12 + mq2_Azero2_coeff * Sqr(Azero);
   ml2 = ml2_m02_coeff * m0Sqr + ml2_m122_coeff * Sqr(m12)
      + ml2_Azerom12_coeff * Azero * m12 + ml2_Azero2_coeff * Sqr(Azero);
   mHd2 = mHd2_m02_coeff * m0Sqr + mHd2_m122_coeff * Sqr(m12)
      + mHd2_Azerom12_coeff * Azero * m12 + mHd2_Azero2_coeff * Sqr(Azero);
   mHu2 = mHu2_m02_coeff * m0Sqr + mHu2_m122_coeff * Sqr(m12)
      + mHu2_Azerom12_coeff * Azero * m12 + mHu2_Azero2_coeff * Sqr(Azero);
   md2 = md2_m02_coeff * m0Sqr + md2_m122_coeff * Sqr(m12)
      + md2_Azerom12_coeff * Azero * m12 + md2_Azero2_coeff * Sqr(Azero);
   mu2 = mu2_m02_coeff * m0Sqr + mu2_m122_coeff * Sqr(m12)
      + mu2_Azerom12_coeff * Azero * m12 + mu2_Azero2_coeff * Sqr(Azero);
   me2 = me2_m02_coeff * m0Sqr + me2_m122_coeff * Sqr(m12)
      + me2_Azerom12_coeff * Azero * m12 + me2_Azero2_coeff * Sqr(Azero);
   MassB = MassB_Azero_coeff * Azero + MassB_m12_coeff * m12;
   MassWB = MassWB_Azero_coeff * Azero + MassWB_m12_coeff * m12;
   MassG = MassG_Azero_coeff * Azero + MassG_m12_coeff * m12;
}

std::ostream& operator<<(std::ostream& ostr, const CMSSM_semianalytic<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
