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

// File generated at Wed 3 Jun 2015 23:53:01

/**
 * @file CSE6SSM_two_scale_model.cpp
 * @brief implementation of the CSE6SSM two scale model class
 *
 * Contains the definition of the CSE6SSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 3 Jun 2015 23:53:01 with FlexibleSUSY
 * 1.1.0 (git commit: v1.1.0) and SARAH 4.5.6 .
 */

#include "CSE6SSM_two_scale_model.hpp"
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

using namespace CSE6SSM_info;

#define CLASSNAME CSE6SSM<Two_scale>

#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

CLASSNAME::CSE6SSM(const CSE6SSM_input_parameters<Two_scale>& input_)
   : Two_scale_model()
   , CSE6SSM_mass_eigenstates()
   , input(input_)
   , number_of_ewsb_iterations(200)
   , ewsb_loop_order(2)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , ewsb_solution(Eigen::Array<double,number_of_tadpole_equations,1>::Zero())
{
}

CLASSNAME::~CSE6SSM()
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

const CSE6SSM_input_parameters<Two_scale>& CLASSNAME::get_input() const
{
   return input;
}

void CLASSNAME::set_input_parameters(const CSE6SSM_input_parameters<Two_scale>& input_)
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
   tadpole[2] = get_ewsb_eq_hh_3();
   tadpole[3] = get_ewsb_eq_hh_4();
   tadpole[4] = get_ewsb_eq_hh_5();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0)) / vd;
      tadpole[1] -= Re(tadpole_hh(1)) / vu;
      tadpole[2] -= Re(tadpole_hh(2)) / vs;
      tadpole[3] -= Re(tadpole_hh(3)) / vsb;
      tadpole[4] -= Re(tadpole_hh(4)) / vphi;

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] -= two_loop_tadpole[0] / vd;
         tadpole[1] -= two_loop_tadpole[1] / vu;
         tadpole[2] -= two_loop_tadpole[2] / vs;

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
   CSE6SSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double s = model->get_input().sInput;

   // N.B. for this version ALambdax is held constant
   double temp = 0.;
   if (is_zero(model->get_TLambdax())) {
      temp = 0.;
   } else if (Abs(model->get_Lambdax()) < std::numeric_limits<double>::epsilon()) {
      throw DivideByZeroError("in CSE6SSM<Two_scale>::tadpole_equations");
   } else {
      temp = model->get_TLambdax() / model->get_Lambdax();
   }
   const double ALambdax = temp;

   // DH:: this updates vsb to correspond to the current
   // estimate for TanTheta, remove it to leave TanTheta
   // free.
   model->set_vs(s * Cos(ArcTan(gsl_vector_get(x, 0))));
   model->set_vsb(s * Sin(ArcTan(gsl_vector_get(x, 0))));
   model->set_Lambdax(gsl_vector_get(x, 1));
   model->set_TLambdax(gsl_vector_get(x, 1) * ALambdax);
   model->set_vphi(gsl_vector_get(x, 2));
   model->set_XiF(gsl_vector_get(x, 3));
   model->set_LXiF(gsl_vector_get(x, 4));


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
   double x_init[number_of_ewsb_equations];
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

   const double s = LOCALINPUT(sInput);

   // N.B. for this version ALambdax is held constant
   double temp = 0.;
   if (is_zero(TLambdax)) {
      temp = 0.;
   } else if (Abs(Lambdax) < std::numeric_limits<double>::epsilon()) {
      throw DivideByZeroError("in CSE6SSM<Two_scale>::solve_ewsb_iteratively_with");
   } else {
      temp = TLambdax / Lambdax;
   }
   const double ALambdax = temp;

   ewsb_solution(0) = solver->get_solution(0);
   ewsb_solution(1) = solver->get_solution(1);
   ewsb_solution(2) = solver->get_solution(2);
   ewsb_solution(3) = solver->get_solution(3);
   ewsb_solution(4) = solver->get_solution(4);

   vs = s * Cos(ArcTan(solver->get_solution(0)));
   vsb = s * Sin(ArcTan(solver->get_solution(0)));
   Lambdax = solver->get_solution(1);
   TLambdax = ALambdax * solver->get_solution(1);
   vphi = solver->get_solution(2);
   XiF = solver->get_solution(3);
   LXiF = solver->get_solution(4);


   return status;
}

int CLASSNAME::check_ewsb_solution(double precision)
{
   double tadpole[number_of_ewsb_equations];

   if (ewsb_loop_order > 0) {
      calculate_DRbar_masses();
   }

   tadpole_equations(tadpole);

   double residual = Abs(tadpole[0]);

   for (std::size_t i = 1; i < number_of_ewsb_equations; ++i) {
      residual += Abs(tadpole[i]);
   } 
   
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

   error = solve_ewsb_iteratively(0);


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

   const auto s = LOCALINPUT(sInput);
   const auto sgnLambdax = Sign(LOCALINPUT(SignLambdax));

   x_init[0] = AbsSqrt((ms2 + 0.0125*Sqr(g1p)*Sqr(QS)*Sqr(s)) 
                       / (msbar2 + 0.0125*Sqr(g1p)*Sqr(QS)*Sqr(s)));

   const double cth = 1. / Sqrt(1. + Sqr(x_init[0]));
   const double sth = cth*x_init[0];
   const double cth2 = Sqr(cth);
   const double sth2 = Sqr(sth);
   const double c2th = cth2 - sth2;

   x_init[1] = sgnLambdax*AbsSqrt(2.0*(mHu2*Sqr(vu) - mHd2*Sqr(vd) + 0.125*Sqr(g2)*Power(vu,4)
                                        + 0.075*Sqr(g1)*Power(vu,4) - 0.125*Sqr(g2)*Power(vd,4)
                                        - 0.075*Sqr(g1)*Power(vd,4) + 0.0125*Sqr(g1p)*
                                        (3.0*Sqr(vd) - 2.0*Sqr(vu))*
                                        (-3.0*Sqr(vd) - 2.0*Sqr(vu) + QS*Sqr(s)*cth2
                                         - QS*Sqr(s)*sth2)) / (Sqr(s)*cth2*(Sqr(vd) - Sqr(vu))));
   
   double lam = x_init[1];

   x_init[2] = (-4. / (s*sth*vu*lam*Conj(Sigmax) + s*sth*vu*Sigmax*Conj(lam)))*
      (mHd2*vd - 0.35355339059327373*s*cth*vu*TLambdax
       - 0.35355339059327373*s*cth*vu*Conj(TLambdax) + 0.5*AbsSqr(lam)*vd*
       (Sqr(vu) + Sqr(s*cth)) + 0.125*Sqr(g2)*Power(vd,3) + 0.075*Sqr(g1)*Power(vd,3)
       - 0.125*Sqr(g2)*vd*Sqr(vu) - 0.075*Sqr(g1)*vd*Sqr(vu) + 0.1125*Sqr(g1p)*
       Power(vd,3) + 0.075*Sqr(g1p)*vd*Sqr(vu) - 0.0375*Sqr(g1p)*QS*vd*Sqr(s)*c2th);

   double phi = x_init[2];

   x_init[3] = ( 2. / (s*cth*Sigmax + s*cth*Conj(Sigmax)))*
      (msbar2*s*sth - 0.35355339059327373*phi*s*cth*MuPhi*Conj(Sigmax)
       - 0.35355339059327373*phi*s*cth*Sigmax*Conj(MuPhi) - 0.35355339059327373*phi*
       s*cth*TSigmax - 0.35355339059327373*phi*s*cth*Conj(TSigmax) + 0.25*phi*vd*vu*
       lam*Conj(Sigmax) + 0.25*phi*vd*vu*Sigmax*Conj(lam) 
       + 0.5*AbsSqr(Sigmax)*s*sth*Sqr(phi) + 0.5*AbsSqr(Sigmax)*s*sth*Sqr(s*cth) - 0.25*s*cth
       *Sqr(phi)*KappaPr*Conj(Sigmax) - 0.25*s*cth*Sqr(phi)*Sigmax*Conj(KappaPr)
       + 0.0375*Sqr(g1p)*QS*s*sth*Sqr(vd) + 0.025*Sqr(g1p)*QS*s*sth*Sqr(vu) - 0.0125*
       Sqr(g1p)*Sqr(QS)*s*sth*Sqr(s*cth) + 0.0125*Sqr(g1p)*Sqr(QS)*Power(s*sth,3));

   double xi = x_init[3];

   x_init[4] = -0.7071067811865475*
      (mphi2*phi + phi*AbsSqr(MuPhi) + Power(phi,3)*AbsSqr(
         KappaPr) + 0.5*phi*BMuPhi + 0.5*phi*Conj(BMuPhi) + 0.7071067811865475*
       MuPhi*Conj(xi) - 0.35355339059327373*MuPhi*Sqr(s)*cth*sth*Conj(Sigmax) 
       - 0.35355339059327373*Sqr(s)*cth*sth*Conj(TSigmax) +
       phi*Conj(xi)*KappaPr - 0.5*phi*Sqr(s)*cth*sth*Conj(Sigmax)*KappaPr + 0.25*vd*s*sth*
       vu*Conj(Sigmax)*lam + 0.7071067811865475*Conj(MuPhi)*xi + phi*Conj(
          KappaPr)*xi - 0.35355339059327373*Sqr(s)*cth*sth*Conj(MuPhi)*Sigmax - 0.5*phi*Sqr(s)*
       cth*sth*Conj(KappaPr)*Sigmax + 0.25*vd*s*sth*vu*Conj(lam)*Sigmax
       + 1.0606601717798212*MuPhi*Conj(KappaPr)*Sqr(phi) +
       0.35355339059327373*Conj(TKappaPr)*Sqr(phi) + 1.0606601717798212*Conj(
          MuPhi)*KappaPr*Sqr(phi) + 0.5*phi*AbsSqr(Sigmax)*Sqr(s*cth) + 0.5*phi*AbsSqr
       (Sigmax)*Sqr(s*sth) + 0.35355339059327373*Sqr(phi)*TKappaPr -
       0.35355339059327373*Sqr(s)*cth*sth*TSigmax);

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
   // save old parameters
   const double vs_old = vs;
   const double vsb_old = vsb;
   const double Lambdax_old = Lambdax;
   const double TLambdax_old = TLambdax;
   // N.B. for this version ALambdax is held constant
   double temp = 0.;
   if (is_zero(TLambdax)) {
      temp = 0.;
   } else if (Abs(Lambdax) < std::numeric_limits<double>::epsilon()) {
      throw DivideByZeroError("in CSE6SSM<Two_scale>::ewsb_step");
   } else {
      temp = TLambdax / Lambdax;
   }
   const double ALambdax_old = temp;
   const double vphi_old = vphi;
   const double XiF_old = XiF;

   int error;
   double tadpole[number_of_ewsb_equations] = { 0. };

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh(0));
      tadpole[1] += Re(tadpole_hh(1));
      tadpole[2] += Re(tadpole_hh(2));
      tadpole[3] += Re(tadpole_hh(3));
      tadpole[4] += Re(tadpole_hh(4));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] += two_loop_tadpole[0];
         tadpole[1] += two_loop_tadpole[1];
         tadpole[2] += two_loop_tadpole[2];

      }
   }

   double TanTheta_new;
   double Lambdax_new;
   double TLambdax_new;
   double vphi_new;
   double XiF_new;
   double LXiF_new;

   const double s = LOCALINPUT(sInput);
   const double SignLambdax = Sign(LOCALINPUT(SignLambdax));

   // update TanTheta
   double delta = 0.5*AbsSqr(Lambdax_old)*Sqr(vs_old)*(Sqr(vd) + Sqr(vu))
      + 0.5*AbsSqr(Sigmax)*Sqr(vs_old)*Sqr(vphi_old) - 0.5*AbsSqr(Sigmax)*Sqr(vsb_old)*Sqr(vphi_old)
      - 0.35355339059327373*vd*vu*vs_old*TLambdax_old - 0.35355339059327373*vd*vu*vs_old*Conj(TLambdax_old) 
      - 0.25*vphi_old*vsb_old*vd*vu*Lambdax_old*Conj(Sigmax) - 0.25*vphi_old*vsb_old*vd*vu*Sigmax*Conj(Lambdax_old)
      - 0.0375*QS*Sqr(g1p)*Sqr(s)*Sqr(vd) - 0.025*QS*Sqr(g1p)*Sqr(s)*Sqr(vu);
   
   if (ewsb_loop_order > 0) {
      delta -= (vs_old * Re(tadpole_hh(2)) - vsb_old * Re(tadpole_hh(3)));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         delta -= vs_old * two_loop_tadpole[2];
      }
   }

   TanTheta_new = AbsSqrt((ms2 + 0.0125*Sqr(g1p)*Sqr(QS)*Sqr(s) + delta/(Sqr(vs)))
                          / (msbar2 + 0.0125*Sqr(g1p)*Sqr(QS)*Sqr(s)));

   const double vs_new = s * Cos(ArcTan(TanTheta_new));
   const double vsb_new = s * Sin(ArcTan(TanTheta_new));

   vs = vs_new;
   vsb = vsb_new;

   double rhs_Lambdax = mHd2*Sqr(vd) - mHu2*Sqr(vu) + 0.125*Sqr(g2)*Power(vd,4)
      + 0.075*Sqr(g1)*Power(vd,4) - 0.125*Sqr(g2)*Power(vu,4) - 0.075*Sqr(g1)*
      Power(vu,4) + 0.1125*Sqr(g1p)*Power(vd,4) - 0.05*Sqr(g1p)*Power(vu,4)
      - 0.0375*QS*Sqr(g1p)*Sqr(vd)*Sqr(vs_new) + 0.0375*QS*Sqr(g1p)*Sqr(vd)*Sqr(vsb_new)
      + 0.025*QS*Sqr(g1p)*Sqr(vu)*Sqr(vs_new) - 0.025*QS*Sqr(g1p)*Sqr(vu)*Sqr(vsb_new);
   
   if (ewsb_loop_order > 0) {
      // DH:: should have error checking here
      rhs_Lambdax -= (vd*Re(tadpole_hh(0)) - vu*Re(tadpole_hh(1)));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         rhs_Lambdax -= (vd*two_loop_tadpole[0] - vu*two_loop_tadpole[1]);
      }
   }

   rhs_Lambdax *= (2. / (Sqr(vs_new)*(Sqr(vu) - Sqr(vd))));

   Lambdax_new = SignLambdax * AbsSqrt(rhs_Lambdax);
   TLambdax_new = Lambdax_new * ALambdax_old;

   Lambdax = Lambdax_new;
   TLambdax = TLambdax_new;

   double rhs_vphi = mHd2*vd - 0.35355339059327373*vs_new*vu*TLambdax_new
      - 0.35355339059327373*vs_new*vu*Conj(TLambdax_new) + 0.5*AbsSqr(Lambdax_new)*vd*
      (Sqr(vu) + Sqr(vs_new)) + 0.125*Sqr(g2)*Power(vd,3) + 0.075*Sqr(g1)*Power(vd,3)
      - 0.125*Sqr(g2)*vd*Sqr(vu) - 0.075*Sqr(g1)*vd*Sqr(vu) + 0.1125*Sqr(g1p)*
      Power(vd,3) + 0.075*Sqr(g1p)*vd*Sqr(vu) - 0.0375*Sqr(g1p)*QS*vd*
      (Sqr(vs_new) - Sqr(vsb_new));

   if (ewsb_loop_order > 0) {
      rhs_vphi -= Re(tadpole_hh(0));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         rhs_vphi -= two_loop_tadpole[0];
      }
   }

   rhs_vphi *= (-4. / (vsb_new*vu*Lambdax_new*Conj(Sigmax) + vsb_new*vu*Sigmax*Conj(Lambdax_new)));

   vphi_new = rhs_vphi;

   vphi = vphi_new;

   double rhs_XiF = msbar2*vsb_new - 0.35355339059327373*vphi_new*vs_new*MuPhi*Conj(Sigmax)
      - 0.35355339059327373*vphi_new*vs_new*Sigmax*Conj(MuPhi) - 0.35355339059327373*vphi_new*
      vs_new*TSigmax - 0.35355339059327373*vphi_new*vs_new*Conj(TSigmax) + 0.25*vphi_new*vd*vu*
      Lambdax_new*Conj(Sigmax) + 0.25*vphi_new*vd*vu*Sigmax*Conj(Lambdax_new) 
      + 0.5*AbsSqr(Sigmax)*vsb_new*Sqr(vphi_new) + 0.5*AbsSqr(Sigmax)*vsb_new*Sqr(vs_new) - 0.25*vs_new
      *Sqr(vphi_new)*KappaPr*Conj(Sigmax) - 0.25*vs_new*Sqr(vphi_new)*Sigmax*Conj(KappaPr)
      + 0.0375*Sqr(g1p)*QS*vsb_new*Sqr(vd) + 0.025*Sqr(g1p)*QS*vsb_new*Sqr(vu) - 0.0125*
      Sqr(g1p)*Sqr(QS)*vsb_new*Sqr(vs_new) + 0.0125*Sqr(g1p)*Sqr(QS)*Power(vsb_new,3);

   if (ewsb_loop_order > 0) {
      rhs_XiF -= Re(tadpole_hh(3));
      if (ewsb_loop_order > 1) {

      }
   }

   rhs_XiF *= ( 2. / (vs_new*Sigmax + vs_new*Conj(Sigmax)));

   XiF_new = rhs_XiF;

   XiF = XiF_new;

   double rhs_LXiF = mphi2*vphi_new + vphi_new*AbsSqr(MuPhi) + Power(vphi_new,3)*AbsSqr(
      KappaPr) + 0.5*vphi_new*BMuPhi + 0.5*vphi_new*Conj(BMuPhi) + 0.7071067811865475*
      MuPhi*Conj(XiF_new) - 0.35355339059327373*MuPhi*vs_new*vsb_new*Conj(Sigmax) 
      - 0.35355339059327373*vs_new*vsb_new*Conj(TSigmax) +
      vphi_new*Conj(XiF_new)*KappaPr - 0.5*vphi_new*vs_new*vsb_new*Conj(Sigmax)*KappaPr + 0.25*vd*vsb_new*
      vu*Conj(Sigmax)*Lambdax_new + 0.7071067811865475*Conj(MuPhi)*XiF_new + vphi_new*Conj(
      KappaPr)*XiF_new - 0.35355339059327373*vs_new*vsb_new*Conj(MuPhi)*Sigmax - 0.5*vphi_new*vs_new*
      vsb_new*Conj(KappaPr)*Sigmax + 0.25*vd*vsb_new*vu*Conj(Lambdax_new)*Sigmax
      + 1.0606601717798212*MuPhi*Conj(KappaPr)*Sqr(vphi_new) +
      0.35355339059327373*Conj(TKappaPr)*Sqr(vphi_new) + 1.0606601717798212*Conj(
      MuPhi)*KappaPr*Sqr(vphi_new) + 0.5*vphi_new*AbsSqr(Sigmax)*Sqr(vs_new) + 0.5*vphi_new*AbsSqr
      (Sigmax)*Sqr(vsb_new) + 0.35355339059327373*Sqr(vphi_new)*TKappaPr -
      0.35355339059327373*vs_new*vsb_new*TSigmax;

   if (ewsb_loop_order > 0) {
      rhs_LXiF -= Re(tadpole_hh(4));
      if (ewsb_loop_order > 1) {

      }
   }

   rhs_LXiF *= -0.7071067811865475;

   LXiF_new = rhs_LXiF;

   // reset old parameters
   vs = vs_old;
   vsb = vsb_old;
   Lambdax = Lambdax_old;
   TLambdax = TLambdax_old;
   vphi = vphi_old;
   XiF = XiF_old;

   const bool isfinite = IsFinite(TanTheta_new) && IsFinite(Lambdax_new) 
      && IsFinite(vphi_new) && IsFinite(XiF_new) && IsFinite(LXiF_new);

   if (isfinite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = TanTheta_new;
      ewsb_parameters[1] = Lambdax_new;
      ewsb_parameters[2] = vphi_new;
      ewsb_parameters[3] = XiF_new;
      ewsb_parameters[4] = LXiF_new;

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
   CSE6SSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double s = INPUT(sInput);

   // N.B. for this version ALambdax is held constant
   double temp = 0.;
   if (is_zero(model->get_TLambdax())) {
      temp = 0.;
   } else if (Abs(model->get_Lambdax()) < std::numeric_limits<double>::epsilon()) {
      throw DivideByZeroError("in CSE6SSM<Two_scale>::ewsb_step");
   } else {
      temp = model->get_TLambdax() / model->get_Lambdax();
   }
   const double ALambdax = temp;

   const double vs = s * Cos(ArcTan(gsl_vector_get(x, 0)));
   const double vsb = s * Sin(ArcTan(gsl_vector_get(x, 0)));
   const double Lambdax = gsl_vector_get(x, 1);
   const double vphi = gsl_vector_get(x, 2);
   const double XiF = gsl_vector_get(x, 3);
   const double LXiF = gsl_vector_get(x, 4);

   model->set_vs(vs);
   model->set_vsb(vsb);
   model->set_Lambdax(Lambdax);
   model->set_TLambdax(Lambdax * ALambdax);
   model->set_vphi(vphi);
   model->set_XiF(XiF);
   model->set_LXiF(LXiF);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { vsb / vs, Lambdax, vphi, XiF, LXiF };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "CSE6SSM (solver type: two_scale)\n"
           "========================================\n";
   CSE6SSM_mass_eigenstates::print(ostr);
}

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
   CSE6SSM_mass_eigenstates::clear();
}

std::string CLASSNAME::name() const
{
   return "CSE6SSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   CSE6SSM_mass_eigenstates::run_to(scale, eps);
}

/**
 * @brief finds the expansion of the low energy soft gaugino masses
 *
 * This function calculates the expansion of the soft gaugino masses
 * in terms of the input parameters m0, m12 and Azero. This is
 * \f$ M_i = p A_0 + q M_{1/2} \f$.
 * 
 * @param soft gaugino mass to calculate for
 * @param low scale to calculate expansion at
 * @param high scale at which input parameters are defined
 * @return array of coefficients (p, q) as defined above 
 */
Eigen::Array<double,2,1> CLASSNAME::get_soft_gaugino_mass_coeffs(CSE6SSM_info::Parameters gaugino_mass,
                                                                 double low_scale, double high_scale) const
{
   const auto m0 = LOCALINPUT(m0);
   const auto m12 = LOCALINPUT(m12);
   const auto Azero = LOCALINPUT(Azero);

   CSE6SSM_soft_parameters run_model;
   run_model.set_loops(get_loops());
   run_model.set_scale(get_scale());
   run_model.set_thresholds(get_thresholds());
   run_model.set(get());

   run_model.run_to(high_scale);

   // generate points
   double m0_pt1 = m0;
   double m12_pt1 = 0.;
   double Azero_pt1 = Azero;

   if (is_zero(Abs(Azero_pt1))) {
      Azero_pt1 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt1, m12_pt1, Azero_pt1);

   run_model.run_to(low_scale);

   double mass_pt1;
   switch(gaugino_mass) {
   case CSE6SSM_info::MassB: {
      mass_pt1 = run_model.get_MassB(); 
      break;
   }
   case CSE6SSM_info::MassWB: {
      mass_pt1 = run_model.get_MassWB();
      break;
   }
   case CSE6SSM_info::MassG: {
      mass_pt1 = run_model.get_MassG();
      break;
   }
   case CSE6SSM_info::MassBp: {
      mass_pt1 = run_model.get_MassBp();
      break;
   }
   default :
      throw UnknownModelParameterError(gaugino_mass);
   }

   run_model.set_scale(get_scale());
   run_model.set(get());

   run_model.run_to(high_scale);

   double m0_pt2 = m0;
   double m12_pt2 = m12;
   double Azero_pt2 = 0.;

   if (is_zero(Abs(m12_pt2))) {
      m12_pt2 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt2, m12_pt2, Azero_pt2);

   run_model.run_to(low_scale);

   double mass_pt2;
   switch(gaugino_mass) {
   case CSE6SSM_info::MassB: {
      mass_pt2 = run_model.get_MassB(); 
      break;
   }
   case CSE6SSM_info::MassWB: {
      mass_pt2 = run_model.get_MassWB();
      break;
   }
   case CSE6SSM_info::MassG: {
      mass_pt2 = run_model.get_MassG();
      break;
   }
   case CSE6SSM_info::MassBp: {
      mass_pt2 = run_model.get_MassBp();
      break;
   }
   default :
      throw UnknownModelParameterError(gaugino_mass);
   }

   const double oneOverdetc = 1.0 / Sqr(Azero_pt1 * m12_pt2 - Azero_pt2 * m12_pt1);
   const double c11 = oneOverdetc * (Sqr(m12_pt1) + Sqr(m12_pt2));
   const double c22 = oneOverdetc * (Sqr(Azero_pt1) + Sqr(Azero_pt2));
   const double c12 = -oneOverdetc * (Azero_pt1 * m12_pt1 + Azero_pt2 * m12_pt2);

   Eigen::Array<double,2,1> coeffs;

   coeffs(0) = c11 * (mass_pt1 * Azero_pt1 + mass_pt2 * Azero_pt2)
      + c12 * (mass_pt1 * m12_pt1 + mass_pt2 * m12_pt2);
   coeffs(1) = c12 * (mass_pt1 * Azero_pt1 + mass_pt2 * Azero_pt2)
      + c22 * (mass_pt1 * m12_pt1 + mass_pt2 * m12_pt2);

   return coeffs;
}

/**
 * @brief finds the expansion of the low energy soft scalar squared masses
 *
 * This function calculates the expansion of the soft scalar squared masses
 * in terms of the input parameters m0, m12 and Azero. This is
 * \f$ m^2 = a m_0^2 + b M_{1/2}^2 + c M_{1/2} A_0 + d A_0^2 \f$.
 * 
 * @param soft scalar mass squared to calculate for
 * @param low scale to calculate expansion at
 * @param high scale at which input parameters are defined
 * @return array of coefficients (a, b, c, d) as defined above 
 */
Eigen::Array<double,4,1> CLASSNAME::get_soft_scalar_mass_coeffs(CSE6SSM_info::Parameters soft_mass,
                                                                double low_scale, double high_scale) const
{
   const auto m0 = LOCALINPUT(m0);
   const auto m12 = LOCALINPUT(m12);
   const auto Azero = LOCALINPUT(Azero);

   CSE6SSM_soft_parameters run_model;
   run_model.set_loops(get_loops());
   run_model.set_scale(get_scale());
   run_model.set_thresholds(get_thresholds());
   run_model.set(get());

   run_model.run_to(high_scale);

   // generate points
   double m0_pt1 = m0;
   double m12_pt1 = 0.;
   double Azero_pt1 = 0.;

   if (is_zero(Abs(m0_pt1))) {
      m0_pt1 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt1, m12_pt1, Azero_pt1);

   run_model.run_to(low_scale);

   double soft_mass_pt1 = get_soft_mass_squared(run_model, soft_mass);

   run_model.set_scale(get_scale());
   run_model.set(get());

   run_model.run_to(high_scale);

   double m0_pt2 = 0.;
   double m12_pt2 = m12;
   double Azero_pt2 = 0.;

   if (is_zero(Abs(m12_pt2))) {
      m12_pt2 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt2, m12_pt2, Azero_pt2);

   run_model.run_to(low_scale);

   double soft_mass_pt2 = get_soft_mass_squared(run_model, soft_mass);

   run_model.set_scale(get_scale());
   run_model.set(get());

   run_model.run_to(high_scale);

   double m0_pt3 = 0.;
   double m12_pt3 = 0.;
   double Azero_pt3 = Azero;

   if (is_zero(Abs(Azero_pt3))) {
      Azero_pt3 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt3, m12_pt3, Azero_pt3);

   run_model.run_to(low_scale);

   double soft_mass_pt3 = get_soft_mass_squared(run_model, soft_mass);

   run_model.set_scale(get_scale());
   run_model.set(get());

   run_model.run_to(high_scale);

   double m0_pt4 = 0.;
   double m12_pt4 = m12;
   double Azero_pt4 = Azero;

   if (is_zero(Abs(m12_pt4))) {
      m12_pt4 = 1.0;
   } 

   if (is_zero(Abs(Azero_pt4))) {
      Azero_pt4 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt4, m12_pt4, Azero_pt4);

   run_model.run_to(low_scale);

   double soft_mass_pt4 = get_soft_mass_squared(run_model, soft_mass);

   Eigen::Matrix<double,4,4> invC;

   invC(0,0) = Power(m0_pt1, 4) + Power(m0_pt2, 4) + Power(m0_pt3, 4) + Power(m0_pt4, 4);
   invC(0,1) = Sqr(m0_pt1) * Sqr(m12_pt1) + Sqr(m0_pt2) * Sqr(m12_pt2) + Sqr(m0_pt3) * Sqr(m12_pt3)
      + Sqr(m0_pt4) * Sqr(m12_pt4);
   invC(0,2) = Sqr(m0_pt1) * m12_pt1 * Azero_pt1 + Sqr(m0_pt2) * m12_pt2 * Azero_pt2
      + Sqr(m0_pt3) * m12_pt3 * Azero_pt3 + Sqr(m0_pt4) * m12_pt4 * Azero_pt4;
   invC(0,3) = Sqr(m0_pt1) * Sqr(Azero_pt1) + Sqr(m0_pt2) * Sqr(Azero_pt2) + Sqr(m0_pt3) * Sqr(Azero_pt3)
      + Sqr(m0_pt4) * Sqr(Azero_pt4);
   invC(1,0) = invC(0,1);
   invC(1,1) = Power(m12_pt1, 4) + Power(m12_pt2, 4) + Power(m12_pt3, 4) + Power(m12_pt4, 4);
   invC(1,2) = Power(m12_pt1, 3) * Azero_pt1 + Power(m12_pt2, 3) * Azero_pt2 + Power(m12_pt3, 3) * Azero_pt3
      + Power(m12_pt4, 3) * Azero_pt4;
   invC(1,3) = Sqr(m12_pt1) * Sqr(Azero_pt1) + Sqr(m12_pt2) * Sqr(Azero_pt2) + Sqr(m12_pt3) * Sqr(Azero_pt3)
      + Sqr(m12_pt4) * Sqr(Azero_pt4);
   invC(2,0) = invC(0,2);
   invC(2,1) = invC(1,2);
   invC(2,2) = invC(1,3);
   invC(2,3) = Power(Azero_pt1, 3) * m12_pt1 + Power(Azero_pt2, 3) * m12_pt2 + Power(Azero_pt3, 3) * m12_pt3
      + Power(Azero_pt4, 3) * m12_pt4;
   invC(3,0) = invC(0,3);
   invC(3,1) = invC(1,3);
   invC(3,2) = invC(2,3);
   invC(3,3) = Power(Azero_pt1, 4) + Power(Azero_pt2, 4) + Power(Azero_pt3, 4) + Power(Azero_pt4, 4);

   Eigen::Matrix<double,4,4> C = invC.inverse();
   
   Eigen::Matrix<double,4,1> b;

   b(0,0) = soft_mass_pt1 * Sqr(m0_pt1) + soft_mass_pt2 * Sqr(m0_pt2) + soft_mass_pt3 * Sqr(m0_pt3)
      + soft_mass_pt4 * Sqr(m0_pt4);
   b(1,0) = soft_mass_pt1 * Sqr(m12_pt1) + soft_mass_pt2 * Sqr(m12_pt2) + soft_mass_pt3 * Sqr(m12_pt3)
      + soft_mass_pt4 * Sqr(m12_pt4);
   b(2,0) = soft_mass_pt1 * m12_pt1 * Azero_pt1 + soft_mass_pt2 * m12_pt2 * Azero_pt2 
      + soft_mass_pt3 * m12_pt3 * Azero_pt3 + soft_mass_pt4 * m12_pt4 * Azero_pt4;
   b(3,0) = soft_mass_pt1 * Sqr(Azero_pt1) + soft_mass_pt2 * Sqr(Azero_pt2) + soft_mass_pt3 * Sqr(Azero_pt3)
      + soft_mass_pt4 * Sqr(Azero_pt4);

   return (C * b).array();
}

double CLASSNAME::get_soft_mass_squared(const CSE6SSM_soft_parameters& model, CSE6SSM_info::Parameters soft_mass) const
{
   double mass2;
   switch(soft_mass) {
   case CSE6SSM_info::mq200: {
      mass2 = model.get_mq2(0,0);
      break;
   }

   case CSE6SSM_info::mq201: {
      mass2 = model.get_mq2(0,1);
      break;
   }

   case CSE6SSM_info::mq202: {
      mass2 = model.get_mq2(0,2);
      break;
   }

   case CSE6SSM_info::mq210: {
      mass2 = model.get_mq2(1,0);
      break;
   }

   case CSE6SSM_info::mq211: {
      mass2 = model.get_mq2(1,1);
      break;
   }

   case CSE6SSM_info::mq212: {
      mass2 = model.get_mq2(1,2);
      break;
   }

   case CSE6SSM_info::mq220: {
      mass2 = model.get_mq2(2,0);
      break;
   }

   case CSE6SSM_info::mq221: {
      mass2 = model.get_mq2(2,1);
      break;
   }

   case CSE6SSM_info::mq222: {
      mass2 = model.get_mq2(2,2);
      break;
   }

   case CSE6SSM_info::ml200: {
      mass2 = model.get_ml2(0,0);
      break;
   }

   case CSE6SSM_info::ml201: {
      mass2 = model.get_ml2(0,1);
      break;
   }

   case CSE6SSM_info::ml202: {
      mass2 = model.get_ml2(0,2);
      break;
   }

   case CSE6SSM_info::ml210: {
      mass2 = model.get_ml2(1,0);
      break;
   }

   case CSE6SSM_info::ml211: {
      mass2 = model.get_ml2(1,1);
      break;
   }

   case CSE6SSM_info::ml212: {
      mass2 = model.get_ml2(1,2);
      break;
   }

   case CSE6SSM_info::ml220: {
      mass2 = model.get_ml2(2,0);
      break;
   }

   case CSE6SSM_info::ml221: {
      mass2 = model.get_ml2(2,1);
      break;
   }

   case CSE6SSM_info::ml222: {
      mass2 = model.get_ml2(2,2);
      break;
   }

   case CSE6SSM_info::md200: {
      mass2 = model.get_md2(0,0);
      break;
   }

   case CSE6SSM_info::md201: {
      mass2 = model.get_md2(0,1);
      break;
   }

   case CSE6SSM_info::md202: {
      mass2 = model.get_md2(0,2);
      break;
   }

   case CSE6SSM_info::md210: {
      mass2 = model.get_md2(1,0);
      break;
   }

   case CSE6SSM_info::md211: {
      mass2 = model.get_md2(1,1);
      break;
   }

   case CSE6SSM_info::md212: {
      mass2 = model.get_md2(1,2);
      break;
   }

   case CSE6SSM_info::md220: {
      mass2 = model.get_md2(2,0);
      break;
   }

   case CSE6SSM_info::md221: {
      mass2 = model.get_md2(2,1);
      break;
   }

   case CSE6SSM_info::md222: {
      mass2 = model.get_md2(2,2);
      break;
   }

   case CSE6SSM_info::mu200: {
      mass2 = model.get_mu2(0,0);
      break;
   }

   case CSE6SSM_info::mu201: {
      mass2 = model.get_mu2(0,1);
      break;
   }

   case CSE6SSM_info::mu202: {
      mass2 = model.get_mu2(0,2);
      break;
   }

   case CSE6SSM_info::mu210: {
      mass2 = model.get_mu2(1,0);
      break;
   }

   case CSE6SSM_info::mu211: {
      mass2 = model.get_mu2(1,1);
      break;
   }

   case CSE6SSM_info::mu212: {
      mass2 = model.get_mu2(1,2);
      break;
   }

   case CSE6SSM_info::mu220: {
      mass2 = model.get_mu2(2,0);
      break;
   }

   case CSE6SSM_info::mu221: {
      mass2 = model.get_mu2(2,1);
      break;
   }

   case CSE6SSM_info::mu222: {
      mass2 = model.get_mu2(2,2);
      break;
   }

   case CSE6SSM_info::me200: {
      mass2 = model.get_me2(0,0);
      break;
   }

   case CSE6SSM_info::me201: {
      mass2 = model.get_me2(0,1);
      break;
   }

   case CSE6SSM_info::me202: {
      mass2 = model.get_me2(0,2);
      break;
   }

   case CSE6SSM_info::me210: {
      mass2 = model.get_me2(1,0);
      break;
   }

   case CSE6SSM_info::me211: {
      mass2 = model.get_me2(1,1);
      break;
   }

   case CSE6SSM_info::me212: {
      mass2 = model.get_me2(1,2);
      break;
   }

   case CSE6SSM_info::me220: {
      mass2 = model.get_me2(2,0);
      break;
   }

   case CSE6SSM_info::me221: {
      mass2 = model.get_me2(2,1);
      break;
   }

   case CSE6SSM_info::me222: {
      mass2 = model.get_me2(2,2);
      break;
   }

   case CSE6SSM_info::mHd2: {
      mass2 = model.get_mHd2();
      break;
   }

   case CSE6SSM_info::mHu2: {
      mass2 = model.get_mHu2();
      break;
   }

   case CSE6SSM_info::ms2: {
      mass2 = model.get_ms2();
      break;
   }

   case CSE6SSM_info::msbar2: {
      mass2 = model.get_msbar2();
      break;
   }

   case CSE6SSM_info::mH1I200: {
      mass2 = model.get_mH1I2(0,0);
      break;
   }

   case CSE6SSM_info::mH1I201: {
      mass2 = model.get_mH1I2(0,1);
      break;
   }

   case CSE6SSM_info::mH1I210: {
      mass2 = model.get_mH1I2(1,0);
      break;
   }

   case CSE6SSM_info::mH1I211: {
      mass2 = model.get_mH1I2(1,1);
      break;
   }

   case CSE6SSM_info::mH2I200: {
      mass2 = model.get_mH2I2(0,0);
      break;
   }

   case CSE6SSM_info::mH2I201: {
      mass2 = model.get_mH2I2(0,1);
      break;
   }

   case CSE6SSM_info::mH2I210: {
      mass2 = model.get_mH2I2(1,0);
      break;
   }

   case CSE6SSM_info::mH2I211: {
      mass2 = model.get_mH2I2(1,1);
      break;
   }

   case CSE6SSM_info::mSI200: {
      mass2 = model.get_mSI2(0,0);
      break;
   }

   case CSE6SSM_info::mSI201: {
      mass2 = model.get_mSI2(0,1);
      break;
   }

   case CSE6SSM_info::mSI202: {
      mass2 = model.get_mSI2(0,2);
      break;
   }

   case CSE6SSM_info::mSI210: {
      mass2 = model.get_mSI2(1,0);
      break;
   }

   case CSE6SSM_info::mSI211: {
      mass2 = model.get_mSI2(1,1);
      break;
   }

   case CSE6SSM_info::mSI212: {
      mass2 = model.get_mSI2(1,2);
      break;
   }

   case CSE6SSM_info::mSI220: {
      mass2 = model.get_mSI2(2,0);
      break;
   }

   case CSE6SSM_info::mSI221: {
      mass2 = model.get_mSI2(2,1);
      break;
   }

   case CSE6SSM_info::mSI222: {
      mass2 = model.get_mSI2(2,2);
      break;
   }

   case CSE6SSM_info::mDx200: {
      mass2 = model.get_mDx2(0,0);
      break;
   }

   case CSE6SSM_info::mDx201: {
      mass2 = model.get_mDx2(0,1);
      break;
   }

   case CSE6SSM_info::mDx202: {
      mass2 = model.get_mDx2(0,2);
      break;
   }

   case CSE6SSM_info::mDx210: {
      mass2 = model.get_mDx2(1,0);
      break;
   }

   case CSE6SSM_info::mDx211: {
      mass2 = model.get_mDx2(1,1);
      break;
   }

   case CSE6SSM_info::mDx212: {
      mass2 = model.get_mDx2(1,2);
      break;
   }

   case CSE6SSM_info::mDx220: {
      mass2 = model.get_mDx2(2,0);
      break;
   }

   case CSE6SSM_info::mDx221: {
      mass2 = model.get_mDx2(2,1);
      break;
   }

   case CSE6SSM_info::mDx222: {
      mass2 = model.get_mDx2(2,2);
      break;
   }

   case CSE6SSM_info::mDxbar200: {
      mass2 = model.get_mDxbar2(0,0);
      break;
   }

   case CSE6SSM_info::mDxbar201: {
      mass2 = model.get_mDxbar2(0,1);
      break;
   }

   case CSE6SSM_info::mDxbar202: {
      mass2 = model.get_mDxbar2(0,2);
      break;
   }

   case CSE6SSM_info::mDxbar210: {
      mass2 = model.get_mDxbar2(1,0);
      break;
   }

   case CSE6SSM_info::mDxbar211: {
      mass2 = model.get_mDxbar2(1,1);
      break;
   }

   case CSE6SSM_info::mDxbar212: {
      mass2 = model.get_mDxbar2(1,2);
      break;
   }

   case CSE6SSM_info::mDxbar220: {
      mass2 = model.get_mDxbar2(2,0);
      break;
   }

   case CSE6SSM_info::mDxbar221: {
      mass2 = model.get_mDxbar2(2,1);
      break;
   }

   case CSE6SSM_info::mDxbar222: {
      mass2 = model.get_mDxbar2(2,2);
      break;
   }

   case CSE6SSM_info::mHp2: {
      mass2 = model.get_mHp2();
      break;
   }

   case CSE6SSM_info::mHpbar2: {
      mass2 = model.get_mHpbar2();
      break;
   }

   case CSE6SSM_info::mphi2: {
      mass2 = model.get_mphi2();
      break;
   }

   default :
      throw UnknownModelParameterError(soft_mass);
   }

   return mass2;
}

/**
 * @brief finds the expansion of the low energy soft trilinears
 *
 * This function calculates the expansion of the soft trilinears
 * in terms of the input parameters m0, m12 and Azero. This is
 * \f$ T_i = e A_0 + f M_{1/2} \f$.
 * 
 * @param soft trilinear to calculate for
 * @param low scale to calculate expansion at
 * @param high scale at which input parameters are defined
 * @return array of coefficients (e, f) as defined above 
 */
Eigen::Array<double,2,1> CLASSNAME::get_soft_trilinear_coeffs(CSE6SSM_info::Parameters soft_trilinear, 
                                                              double low_scale, double high_scale) const
{
   const auto m0 = LOCALINPUT(m0);
   const auto m12 = LOCALINPUT(m12);
   const auto Azero = LOCALINPUT(Azero);

   CSE6SSM_soft_parameters run_model;
   run_model.set_loops(get_loops());
   run_model.set_scale(get_scale());
   run_model.set_thresholds(get_thresholds());
   run_model.set(get());

   run_model.run_to(high_scale);

   // generate points
   double m0_pt1 = m0;
   double m12_pt1 = 0.;
   double Azero_pt1 = Azero;

   if (is_zero(Abs(Azero_pt1))) {
      Azero_pt1 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt1, m12_pt1, Azero_pt1);

   run_model.run_to(low_scale);

   double trilinear_pt1 = get_soft_trilinear(run_model, soft_trilinear);

   run_model.set_scale(get_scale());
   run_model.set(get());

   run_model.run_to(high_scale);

   double m0_pt2 = m0;
   double m12_pt2 = m12;
   double Azero_pt2 = 0.;

   if (is_zero(Abs(m12_pt2))) {
      m12_pt2 = 1.0;
   } 

   set_pars_at_high_scale(run_model, m0_pt2, m12_pt2, Azero_pt2);

   run_model.run_to(low_scale);

   double trilinear_pt2 = get_soft_trilinear(run_model, soft_trilinear);

   const double oneOverdetc = 1.0 / Sqr(Azero_pt1 * m12_pt2 - Azero_pt2 * m12_pt1);
   const double c11 = oneOverdetc * (Sqr(m12_pt1) + Sqr(m12_pt2));
   const double c22 = oneOverdetc * (Sqr(Azero_pt1) + Sqr(Azero_pt2));
   const double c12 = -oneOverdetc * (Azero_pt1 * m12_pt1 + Azero_pt2 * m12_pt2);

   Eigen::Array<double,2,1> coeffs;

   coeffs(0) = c11 * (trilinear_pt1 * Azero_pt1 + trilinear_pt2 * Azero_pt2)
      + c12 * (trilinear_pt1 * m12_pt1 + trilinear_pt2 * m12_pt2);
   coeffs(1) = c12 * (trilinear_pt1 * Azero_pt1 + trilinear_pt2 * Azero_pt2)
      + c22 * (trilinear_pt1 * m12_pt1 + trilinear_pt2 * m12_pt2);

   return coeffs;
}

double CLASSNAME::get_soft_trilinear(const CSE6SSM_soft_parameters& model, CSE6SSM_info::Parameters soft_trilinear) const
{
   double trilinear = 0.;

   switch(soft_trilinear) {
   case CSE6SSM_info::TYe00: {
      trilinear = model.get_TYe(0,0);
      break;
   }

   case CSE6SSM_info::TYe01: {
      trilinear = model.get_TYe(0,1);
      break;
   }

   case CSE6SSM_info::TYe02: {
      trilinear = model.get_TYe(0,2);
      break;
   }

   case CSE6SSM_info::TYe10: {
      trilinear = model.get_TYe(1,0);
      break;
   }

   case CSE6SSM_info::TYe11: {
      trilinear = model.get_TYe(1,1);
      break;
   }

   case CSE6SSM_info::TYe12: {
      trilinear = model.get_TYe(1,2);
      break;
   }

   case CSE6SSM_info::TYe20: {
      trilinear = model.get_TYe(2,0);
      break;
   }

   case CSE6SSM_info::TYe21: {
      trilinear = model.get_TYe(2,1);
      break;
   }

   case CSE6SSM_info::TYe22: {
      trilinear = model.get_TYe(2,2);
      break;
   }

   case CSE6SSM_info::TYd00: {
      trilinear = model.get_TYd(0,0);
      break;
   }

   case CSE6SSM_info::TYd01: {
      trilinear = model.get_TYd(0,1);
      break;
   }

   case CSE6SSM_info::TYd02: {
      trilinear = model.get_TYd(0,2);
      break;
   }

   case CSE6SSM_info::TYd10: {
      trilinear = model.get_TYd(1,0);
      break;
   }

   case CSE6SSM_info::TYd11: {
      trilinear = model.get_TYd(1,1);
      break;
   }

   case CSE6SSM_info::TYd12: {
      trilinear = model.get_TYd(1,2);
      break;
   }

   case CSE6SSM_info::TYd20: {
      trilinear = model.get_TYd(2,0);
      break;
   }

   case CSE6SSM_info::TYd21: {
      trilinear = model.get_TYd(2,1);
      break;
   }

   case CSE6SSM_info::TYd22: {
      trilinear = model.get_TYd(2,2);
      break;
   }

   case CSE6SSM_info::TYu00: {
      trilinear = model.get_TYu(0,0);
      break;
   }

   case CSE6SSM_info::TYu01: {
      trilinear = model.get_TYu(0,1);
      break;
   }

   case CSE6SSM_info::TYu02: {
      trilinear = model.get_TYu(0,2);
      break;
   }

   case CSE6SSM_info::TYu10: {
      trilinear = model.get_TYu(1,0);
      break;
   }

   case CSE6SSM_info::TYu11: {
      trilinear = model.get_TYu(1,1);
      break;
   }

   case CSE6SSM_info::TYu12: {
      trilinear = model.get_TYu(1,2);
      break;
   }

   case CSE6SSM_info::TYu20: {
      trilinear = model.get_TYu(2,0);
      break;
   }

   case CSE6SSM_info::TYu21: {
      trilinear = model.get_TYu(2,1);
      break;
   }

   case CSE6SSM_info::TYu22: {
      trilinear = model.get_TYu(2,2);
      break;
   }

   case CSE6SSM_info::TKappa00: {
      trilinear = model.get_TKappa(0,0);
      break;
   }

   case CSE6SSM_info::TKappa01: {
      trilinear = model.get_TKappa(0,1);
      break;
   }

   case CSE6SSM_info::TKappa02: {
      trilinear = model.get_TKappa(0,2);
      break;
   }

   case CSE6SSM_info::TKappa10: {
      trilinear = model.get_TKappa(1,0);
      break;
   }

   case CSE6SSM_info::TKappa11: {
      trilinear = model.get_TKappa(1,1);
      break;
   }

   case CSE6SSM_info::TKappa12: {
      trilinear = model.get_TKappa(1,2);
      break;
   }

   case CSE6SSM_info::TKappa20: {
      trilinear = model.get_TKappa(2,0);
      break;
   }

   case CSE6SSM_info::TKappa21: {
      trilinear = model.get_TKappa(2,1);
      break;
   }

   case CSE6SSM_info::TKappa22: {
      trilinear = model.get_TKappa(2,2);
      break;
   }

   case CSE6SSM_info::TLambda1200: {
      trilinear = model.get_TLambda12(0,0);
      break;
   }

   case CSE6SSM_info::TLambda1201: {
      trilinear = model.get_TLambda12(0,1);
      break;
   }

   case CSE6SSM_info::TLambda1210: {
      trilinear = model.get_TLambda12(1,0);
      break;
   }

   case CSE6SSM_info::TLambda1211: {
      trilinear = model.get_TLambda12(1,1);
      break;
   }

   case CSE6SSM_info::TLambdax: {
      trilinear = model.get_TLambdax();
      break;
   }

   case CSE6SSM_info::ThE00: {
      trilinear = model.get_ThE(0,0);
      break;
   }

   case CSE6SSM_info::ThE01: {
      trilinear = model.get_ThE(0,1);
      break;
   }

   case CSE6SSM_info::ThE10: {
      trilinear = model.get_ThE(1,0);
      break;
   }

   case CSE6SSM_info::ThE11: {
      trilinear = model.get_ThE(1,1);
      break;
   }

   case CSE6SSM_info::ThE20: {
      trilinear = model.get_ThE(2,0);
      break;
   }

   case CSE6SSM_info::ThE21: {
      trilinear = model.get_ThE(2,1);
      break;
   }

   case CSE6SSM_info::TSigmaL: {
      trilinear = model.get_TSigmaL();
      break;
   }

   case CSE6SSM_info::TKappaPr: {
      trilinear = model.get_TKappaPr();
      break;
   }

   case CSE6SSM_info::TSigmax: {
      trilinear = model.get_TSigmax();
      break;
   }

   case CSE6SSM_info::TgD00: {
      trilinear = model.get_TgD(0,0);
      break;
   }

   case CSE6SSM_info::TgD01: {
      trilinear = model.get_TgD(0,1);
      break;
   }

   case CSE6SSM_info::TgD02: {
      trilinear = model.get_TgD(0,2);
      break;
   }

   case CSE6SSM_info::TgD10: {
      trilinear = model.get_TgD(1,0);
      break;
   }

   case CSE6SSM_info::TgD11: {
      trilinear = model.get_TgD(1,1);
      break;
   }

   case CSE6SSM_info::TgD12: {
      trilinear = model.get_TgD(1,2);
      break;
   }

   case CSE6SSM_info::TgD20: {
      trilinear = model.get_TgD(2,0);
      break;
   }

   case CSE6SSM_info::TgD21: {
      trilinear = model.get_TgD(2,1);
      break;
   }

   case CSE6SSM_info::TgD22: {
      trilinear = model.get_TgD(2,2);
      break;
   }

   case CSE6SSM_info::Tfu00: {
      trilinear = model.get_Tfu(0,0);
      break;
   }

   case CSE6SSM_info::Tfu01: {
      trilinear = model.get_Tfu(0,1);
      break;
   }

   case CSE6SSM_info::Tfu10: {
      trilinear = model.get_Tfu(1,0);
      break;
   }

   case CSE6SSM_info::Tfu11: {
      trilinear = model.get_Tfu(1,1);
      break;
   }

   case CSE6SSM_info::Tfu20: {
      trilinear = model.get_Tfu(2,0);
      break;
   }

   case CSE6SSM_info::Tfu21: {
      trilinear = model.get_Tfu(2,1);
      break;
   }

   case CSE6SSM_info::Tfd00: {
      trilinear = model.get_Tfd(0,0);
      break;
   }

   case CSE6SSM_info::Tfd01: {
      trilinear = model.get_Tfd(0,1);
      break;
   }

   case CSE6SSM_info::Tfd10: {
      trilinear = model.get_Tfd(1,0);
      break;
   }

   case CSE6SSM_info::Tfd11: {
      trilinear = model.get_Tfd(1,1);
      break;
   }

   case CSE6SSM_info::Tfd20: {
      trilinear = model.get_Tfd(2,0);
      break;
   }

   case CSE6SSM_info::Tfd21: {
      trilinear = model.get_Tfd(2,1);
      break;
   }

   default :
      throw UnknownModelParameterError(soft_trilinear);
   }

   return trilinear;
}

void CLASSNAME::set_pars_at_high_scale(CSE6SSM_soft_parameters & model, double m0, double m12, double Azero) const
{
   model.set_TYe(Azero*model.get_Ye());
   model.set_TYd(Azero*model.get_Yd());
   model.set_TYu(Azero*model.get_Yu());
   model.set_TKappaPr(Azero*model.get_KappaPr());
   model.set_TSigmax(Azero*model.get_Sigmax());
   model.set_ThE(Azero*model.get_hE());
   model.set_TSigmaL(Azero*model.get_SigmaL());
   model.set_TgD(Azero*model.get_gD());
   model.set_Tfu(Azero*model.get_fu());
   model.set_Tfd(Azero*model.get_fd());
   model.set_TKappa(Azero*model.get_Kappa());
   model.set_TLambda12(Azero*model.get_Lambda12());
   model.set_TLambdax(Azero*model.get_Lambdax());
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_ms2(Sqr(m0));
   model.set_msbar2(Sqr(m0));
   model.set_mphi2(Sqr(m0));
   model.set_mHp2(Sqr(m0));
   model.set_mHpbar2(Sqr(m0));
   model.set_mH1I2(Sqr(m0)*UNITMATRIX(2));
   model.set_mH2I2(Sqr(m0)*UNITMATRIX(2));
   model.set_mSI2(Sqr(m0)*UNITMATRIX(3));
   model.set_mq2(Sqr(m0)*UNITMATRIX(3));
   model.set_ml2(Sqr(m0)*UNITMATRIX(3));
   model.set_md2(Sqr(m0)*UNITMATRIX(3));
   model.set_mu2(Sqr(m0)*UNITMATRIX(3));
   model.set_me2(Sqr(m0)*UNITMATRIX(3));
   model.set_mDx2(Sqr(m0)*UNITMATRIX(3));
   model.set_mDxbar2(Sqr(m0)*UNITMATRIX(3));
   model.set_MassB(m12);
   model.set_MassWB(m12);
   model.set_MassG(m12);
   model.set_MassBp(m12);
}

double CLASSNAME::get_parameter(unsigned parameter) const
{
   if (parameter >= CSE6SSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CSE6SSM_info::Yd00:
      return Yd(0,0);
   case CSE6SSM_info::Yd01:
      return Yd(0,1);
   case CSE6SSM_info::Yd02:
      return Yd(0,2);
   case CSE6SSM_info::Yd10:
      return Yd(1,0); 
   case CSE6SSM_info::Yd11:
      return Yd(1,1); 
   case CSE6SSM_info::Yd12:
      return Yd(1,2);
   case CSE6SSM_info::Yd20:
      return Yd(2,0);
   case CSE6SSM_info::Yd21:
      return Yd(2,1);
   case CSE6SSM_info::Yd22:
      return Yd(2,2);
   case CSE6SSM_info::hE00:
      return hE(0,0);
   case CSE6SSM_info::hE01:
      return hE(0,1);
   case CSE6SSM_info::hE10:
      return hE(1,0); 
   case CSE6SSM_info::hE11:
      return hE(1,1); 
   case CSE6SSM_info::hE20:
      return hE(2,0); 
   case CSE6SSM_info::hE21:
      return hE(2,1);
   case CSE6SSM_info::Ye00:
      return Ye(0,0); 
   case CSE6SSM_info::Ye01:
      return Ye(0,1);
   case CSE6SSM_info::Ye02:
      return Ye(0,2); 
   case CSE6SSM_info::Ye10:
      return Ye(1,0); 
   case CSE6SSM_info::Ye11:
      return Ye(1,1); 
   case CSE6SSM_info::Ye12:
      return Ye(1,2);
   case CSE6SSM_info::Ye20:
      return Ye(2,0); 
   case CSE6SSM_info::Ye21:
      return Ye(2,1); 
   case CSE6SSM_info::Ye22:
      return Ye(2,2); 
   case CSE6SSM_info::SigmaL:
      return SigmaL; 
   case CSE6SSM_info::KappaPr:
      return KappaPr; 
   case CSE6SSM_info::Sigmax:
      return Sigmax; 
   case CSE6SSM_info::gD00:
      return gD(0,0); 
   case CSE6SSM_info::gD01:
      return gD(0,1); 
   case CSE6SSM_info::gD02:
      return gD(0,2); 
   case CSE6SSM_info::gD10:
      return gD(1,0); 
   case CSE6SSM_info::gD11:
      return gD(1,1);
   case CSE6SSM_info::gD12:
      return gD(1,2); 
   case CSE6SSM_info::gD20:
      return gD(2,0);
   case CSE6SSM_info::gD21:
      return gD(2,1); 
   case CSE6SSM_info::gD22:
      return gD(2,2); 
   case CSE6SSM_info::Kappa00:
      return Kappa(0,0); 
   case CSE6SSM_info::Kappa01:
      return Kappa(0,1); 
   case CSE6SSM_info::Kappa02:
      return Kappa(0,2); 
   case CSE6SSM_info::Kappa10:
      return Kappa(1,0); 
   case CSE6SSM_info::Kappa11:
      return Kappa(1,1); 
   case CSE6SSM_info::Kappa12:
      return Kappa(1,2);
   case CSE6SSM_info::Kappa20:
      return Kappa(2,0); 
   case CSE6SSM_info::Kappa21:
      return Kappa(2,1); 
   case CSE6SSM_info::Kappa22:
      return Kappa(2,2); 
   case CSE6SSM_info::Lambda1200:
      return Lambda12(0,0); 
   case CSE6SSM_info::Lambda1201:
      return Lambda12(0,1); 
   case CSE6SSM_info::Lambda1210:
      return Lambda12(1,0); 
   case CSE6SSM_info::Lambda1211:
      return Lambda12(1,1);
   case CSE6SSM_info::Lambdax:
      return Lambdax;
   case CSE6SSM_info::fu00:
      return fu(0,0);
   case CSE6SSM_info::fu01:
      return fu(0,1); 
   case CSE6SSM_info::fu10:
      return fu(1,0);
   case CSE6SSM_info::fu11:
      return fu(1,1); 
   case CSE6SSM_info::fu20:
      return fu(2,0); 
   case CSE6SSM_info::fu21:
      return fu(2,1); 
   case CSE6SSM_info::fd00:
      return fd(0,0); 
   case CSE6SSM_info::fd01:
      return fd(0,1);
   case CSE6SSM_info::fd10:
      return fd(1,0); 
   case CSE6SSM_info::fd11:
      return fd(1,1); 
   case CSE6SSM_info::fd20:
      return fd(2,0);
   case CSE6SSM_info::fd21:
      return fd(2,1); 
   case CSE6SSM_info::Yu00:
      return Yu(0,0); 
   case CSE6SSM_info::Yu01:
      return Yu(0,1); 
   case CSE6SSM_info::Yu02:
      return Yu(0,2); 
   case CSE6SSM_info::Yu10:
      return Yu(1,0); 
   case CSE6SSM_info::Yu11:
      return Yu(1,1); 
   case CSE6SSM_info::Yu12:
      return Yu(1,2); 
   case CSE6SSM_info::Yu20:
      return Yu(2,0); 
   case CSE6SSM_info::Yu21:
      return Yu(2,1); 
   case CSE6SSM_info::Yu22:
      return Yu(2,2); 
   case CSE6SSM_info::MuPr:
      return MuPr;
   case CSE6SSM_info::MuPhi:
      return MuPhi; 
   case CSE6SSM_info::XiF:
      return XiF;
   case CSE6SSM_info::g1:
      return g1; 
   case CSE6SSM_info::g2:
      return g2;
   case CSE6SSM_info::g3:
      return g3; 
   case CSE6SSM_info::g1p:
      return g1p; 
   case CSE6SSM_info::vd:
      return vd; 
   case CSE6SSM_info::vu:
      return vu; 
   case CSE6SSM_info::vs:
      return vs; 
   case CSE6SSM_info::vsb:
      return vsb; 
   case CSE6SSM_info::vphi:
      return vphi; 
   case CSE6SSM_info::QS:
      return QS;
   case CSE6SSM_info::TYd00:
      return TYd(0,0); 
   case CSE6SSM_info::TYd01:
      return TYd(0,1); 
   case CSE6SSM_info::TYd02:
      return TYd(0,2); 
   case CSE6SSM_info::TYd10:
      return TYd(1,0); 
   case CSE6SSM_info::TYd11:
      return TYd(1,1);
   case CSE6SSM_info::TYd12:
      return TYd(1,2); 
   case CSE6SSM_info::TYd20:
      return TYd(2,0); 
   case CSE6SSM_info::TYd21:
      return TYd(2,1);
   case CSE6SSM_info::TYd22:
      return TYd(2,2); 
   case CSE6SSM_info::ThE00:
      return ThE(0,0); 
   case CSE6SSM_info::ThE01:
      return ThE(0,1); 
   case CSE6SSM_info::ThE10:
      return ThE(1,0); 
   case CSE6SSM_info::ThE11:
      return ThE(1,1); 
   case CSE6SSM_info::ThE20:
      return ThE(2,0); 
   case CSE6SSM_info::ThE21:
      return ThE(2,1); 
   case CSE6SSM_info::TYe00:
      return TYe(0,0);
   case CSE6SSM_info::TYe01:
      return TYe(0,1); 
   case CSE6SSM_info::TYe02:
      return TYe(0,2); 
   case CSE6SSM_info::TYe10:
      return TYe(1,0); 
   case CSE6SSM_info::TYe11:
      return TYe(1,1); 
   case CSE6SSM_info::TYe12:
      return TYe(1,2); 
   case CSE6SSM_info::TYe20:
      return TYe(2,0); 
   case CSE6SSM_info::TYe21:
      return TYe(2,1); 
   case CSE6SSM_info::TYe22:
      return TYe(2,2); 
   case CSE6SSM_info::TSigmaL:
      return TSigmaL; 
   case CSE6SSM_info::TKappaPr:
      return TKappaPr;
   case CSE6SSM_info::TSigmax:
      return TSigmax;
   case CSE6SSM_info::TgD00:
      return TgD(0,0); 
   case CSE6SSM_info::TgD01:
      return TgD(0,1); 
   case CSE6SSM_info::TgD02:
      return TgD(0,2); 
   case CSE6SSM_info::TgD10:
      return TgD(1,0); 
   case CSE6SSM_info::TgD11:
      return TgD(1,1); 
   case CSE6SSM_info::TgD12:
      return TgD(1,2); 
   case CSE6SSM_info::TgD20:
      return TgD(2,0); 
   case CSE6SSM_info::TgD21:
      return TgD(2,1); 
   case CSE6SSM_info::TgD22:
      return TgD(2,2);
   case CSE6SSM_info::TKappa00:
      return TKappa(0,0); 
   case CSE6SSM_info::TKappa01:
      return TKappa(0,1); 
   case CSE6SSM_info::TKappa02:
      return TKappa(0,2); 
   case CSE6SSM_info::TKappa10:
      return TKappa(1,0); 
   case CSE6SSM_info::TKappa11:
      return TKappa(1,1); 
   case CSE6SSM_info::TKappa12:
      return TKappa(1,2); 
   case CSE6SSM_info::TKappa20:
      return TKappa(2,0);
   case CSE6SSM_info::TKappa21:
      return TKappa(2,1); 
   case CSE6SSM_info::TKappa22:
      return TKappa(2,2);
   case CSE6SSM_info::TLambda1200:
      return TLambda12(0,0); 
   case CSE6SSM_info::TLambda1201:
      return TLambda12(0,1); 
   case CSE6SSM_info::TLambda1210:
      return TLambda12(1,0); 
   case CSE6SSM_info::TLambda1211:
      return TLambda12(1,1);
   case CSE6SSM_info::TLambdax:
      return TLambdax;
   case CSE6SSM_info::Tfu00:
      return Tfu(0,0); 
   case CSE6SSM_info::Tfu01:
      return Tfu(0,1); 
   case CSE6SSM_info::Tfu10:
      return Tfu(1,0);
   case CSE6SSM_info::Tfu11:
      return Tfu(1,1); 
   case CSE6SSM_info::Tfu20:
      return Tfu(2,0); 
   case CSE6SSM_info::Tfu21:
      return Tfu(2,1); 
   case CSE6SSM_info::Tfd00:
      return Tfd(0,0); 
   case CSE6SSM_info::Tfd01:
      return Tfd(0,1); 
   case CSE6SSM_info::Tfd10:
      return Tfd(1,0);
   case CSE6SSM_info::Tfd11:
      return Tfd(1,1); 
   case CSE6SSM_info::Tfd20:
      return Tfd(2,0); 
   case CSE6SSM_info::Tfd21:
      return Tfd(2,1);
   case CSE6SSM_info::TYu00:
      return TYu(0,0); 
   case CSE6SSM_info::TYu01:
      return TYu(0,1); 
   case CSE6SSM_info::TYu02:
      return TYu(0,2); 
   case CSE6SSM_info::TYu10:
      return TYu(1,0); 
   case CSE6SSM_info::TYu11:
      return TYu(1,1); 
   case CSE6SSM_info::TYu12:
      return TYu(1,2); 
   case CSE6SSM_info::TYu20:
      return TYu(2,0); 
   case CSE6SSM_info::TYu21:
      return TYu(2,1);
   case CSE6SSM_info::TYu22:
      return TYu(2,2); 
   case CSE6SSM_info::BMuPr:
      return BMuPr;
   case CSE6SSM_info::BMuPhi:
      return BMuPhi; 
   case CSE6SSM_info::LXiF:
      return LXiF; 
   case CSE6SSM_info::mq200:
      return mq2(0,0); 
   case CSE6SSM_info::mq201:
      return mq2(0,1); 
   case CSE6SSM_info::mq202:
      return mq2(0,2); 
   case CSE6SSM_info::mq210:
      return mq2(1,0); 
   case CSE6SSM_info::mq211:
      return mq2(1,1); 
   case CSE6SSM_info::mq212:
      return mq2(1,2); 
   case CSE6SSM_info::mq220:
      return mq2(2,0);
   case CSE6SSM_info::mq221:
      return mq2(2,1); 
   case CSE6SSM_info::mq222:
      return mq2(2,2); 
   case CSE6SSM_info::ml200:
      return ml2(0,0);
   case CSE6SSM_info::ml201:
      return ml2(0,1); 
   case CSE6SSM_info::ml202:
      return ml2(0,2); 
   case CSE6SSM_info::ml210:
      return ml2(1,0); 
   case CSE6SSM_info::ml211:
      return ml2(1,1); 
   case CSE6SSM_info::ml212:
      return ml2(1,2); 
   case CSE6SSM_info::ml220:
      return ml2(2,0); 
   case CSE6SSM_info::ml221:
      return ml2(2,1);
   case CSE6SSM_info::ml222:
      return ml2(2,2);
   case CSE6SSM_info::mHd2:
      return mHd2; 
   case CSE6SSM_info::mHu2:
      return mHu2; 
   case CSE6SSM_info::md200:
      return md2(0,0); 
   case CSE6SSM_info::md201:
      return md2(0,1); 
   case CSE6SSM_info::md202:
      return md2(0,2); 
   case CSE6SSM_info::md210:
      return md2(1,0); 
   case CSE6SSM_info::md211:
      return md2(1,1); 
   case CSE6SSM_info::md212:
      return md2(1,2);
   case CSE6SSM_info::md220:
      return md2(2,0); 
   case CSE6SSM_info::md221:
      return md2(2,1);
   case CSE6SSM_info::md222:
      return md2(2,2); 
   case CSE6SSM_info::mu200:
      return mu2(0,0); 
   case CSE6SSM_info::mu201:
      return mu2(0,1); 
   case CSE6SSM_info::mu202:
      return mu2(0,2); 
   case CSE6SSM_info::mu210:
      return mu2(1,0); 
   case CSE6SSM_info::mu211:
      return mu2(1,1); 
   case CSE6SSM_info::mu212:
      return mu2(1,2); 
   case CSE6SSM_info::mu220:
      return mu2(2,0); 
   case CSE6SSM_info::mu221:
      return mu2(2,1); 
   case CSE6SSM_info::mu222:
      return mu2(2,2); 
   case CSE6SSM_info::me200:
      return me2(0,0);
   case CSE6SSM_info::me201:
      return me2(0,1); 
   case CSE6SSM_info::me202:
      return me2(0,2); 
   case CSE6SSM_info::me210:
      return me2(1,0); 
   case CSE6SSM_info::me211:
      return me2(1,1); 
   case CSE6SSM_info::me212:
      return me2(1,2); 
   case CSE6SSM_info::me220:
      return me2(2,0);
   case CSE6SSM_info::me221:
      return me2(2,1); 
   case CSE6SSM_info::me222:
      return me2(2,2); 
   case CSE6SSM_info::ms2:
      return ms2; 
   case CSE6SSM_info::msbar2:
      return msbar2;
   case CSE6SSM_info::mH1I200:
      return mH1I2(0,0); 
   case CSE6SSM_info::mH1I201:
      return mH1I2(0,1); 
   case CSE6SSM_info::mH1I210:
      return mH1I2(1,0); 
   case CSE6SSM_info::mH1I211:
      return mH1I2(1,1); 
   case CSE6SSM_info::mH2I200:
      return mH2I2(0,0); 
   case CSE6SSM_info::mH2I201:
      return mH2I2(0,1); 
   case CSE6SSM_info::mH2I210:
      return mH2I2(1,0); 
   case CSE6SSM_info::mH2I211:
      return mH2I2(1,1);
   case CSE6SSM_info::mSI200:
      return mSI2(0,0); 
   case CSE6SSM_info::mSI201:
      return mSI2(0,1); 
   case CSE6SSM_info::mSI202:
      return mSI2(0,2); 
   case CSE6SSM_info::mSI210:
      return mSI2(1,0); 
   case CSE6SSM_info::mSI211:
      return mSI2(1,1); 
   case CSE6SSM_info::mSI212:
      return mSI2(1,2); 
   case CSE6SSM_info::mSI220:
      return mSI2(2,0);
   case CSE6SSM_info::mSI221:
      return mSI2(2,1); 
   case CSE6SSM_info::mSI222:
      return mSI2(2,2);
   case CSE6SSM_info::mDx200:
      return mDx2(0,0); 
   case CSE6SSM_info::mDx201:
      return mDx2(0,1); 
   case CSE6SSM_info::mDx202:
      return mDx2(0,2); 
   case CSE6SSM_info::mDx210:
      return mDx2(1,0); 
   case CSE6SSM_info::mDx211:
      return mDx2(1,1); 
   case CSE6SSM_info::mDx212:
      return mDx2(1,2); 
   case CSE6SSM_info::mDx220:
      return mDx2(2,0); 
   case CSE6SSM_info::mDx221:
      return mDx2(2,1); 
   case CSE6SSM_info::mDx222:
      return mDx2(2,2);
   case CSE6SSM_info::mDxbar200:
      return mDxbar2(0,0); 
   case CSE6SSM_info::mDxbar201:
      return mDxbar2(0,1); 
   case CSE6SSM_info::mDxbar202:
      return mDxbar2(0,2); 
   case CSE6SSM_info::mDxbar210:
      return mDxbar2(1,0); 
   case CSE6SSM_info::mDxbar211:
      return mDxbar2(1,1); 
   case CSE6SSM_info::mDxbar212:
      return mDxbar2(1,2); 
   case CSE6SSM_info::mDxbar220:
      return mDxbar2(2,0);
   case CSE6SSM_info::mDxbar221:
      return mDxbar2(2,1);
   case CSE6SSM_info::mDxbar222:
      return mDxbar2(2,2); 
   case CSE6SSM_info::mHp2:
      return mHp2; 
   case CSE6SSM_info::mHpbar2:
      return mHpbar2; 
   case CSE6SSM_info::mphi2:
      return mphi2; 
   case CSE6SSM_info::MassB:
      return MassB; 
   case CSE6SSM_info::MassWB:
      return MassWB; 
   case CSE6SSM_info::MassG:
      return MassG; 
   case CSE6SSM_info::MassBp:
      return MassBp;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

void CLASSNAME::set_parameter(unsigned parameter, double x)
{
   if (parameter >= CSE6SSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CSE6SSM_info::Yd00:
      Yd(0,0) = x;
      break;
   case CSE6SSM_info::Yd01:
      Yd(0,1) = x;
      break;
   case CSE6SSM_info::Yd02:
      Yd(0,2) = x;
      break;
   case CSE6SSM_info::Yd10:
      Yd(1,0) = x;
      break;
   case CSE6SSM_info::Yd11:
      Yd(1,1) = x;
      break;
   case CSE6SSM_info::Yd12:
      Yd(1,2) = x;
      break;
   case CSE6SSM_info::Yd20:
      Yd(2,0) = x;
      break;
   case CSE6SSM_info::Yd21:
      Yd(2,1) = x;
      break;
   case CSE6SSM_info::Yd22:
      Yd(2,2) = x;
      break;
   case CSE6SSM_info::hE00:
      hE(0,0) = x;
      break;
   case CSE6SSM_info::hE01:
      hE(0,1) = x;
      break;
   case CSE6SSM_info::hE10:
      hE(1,0) = x;
      break;
   case CSE6SSM_info::hE11:
      hE(1,1) = x;
      break;
   case CSE6SSM_info::hE20:
      hE(2,0) = x;
      break;
   case CSE6SSM_info::hE21:
      hE(2,1) = x;
      break;
   case CSE6SSM_info::Ye00:
      Ye(0,0) = x; 
      break;
   case CSE6SSM_info::Ye01:
      Ye(0,1) = x;
      break;
   case CSE6SSM_info::Ye02:
      Ye(0,2) = x;
      break;
   case CSE6SSM_info::Ye10:
      Ye(1,0) = x;
      break;
   case CSE6SSM_info::Ye11:
      Ye(1,1) = x;
      break;
   case CSE6SSM_info::Ye12:
      Ye(1,2) = x;
      break;
   case CSE6SSM_info::Ye20:
      Ye(2,0) = x;
      break;
   case CSE6SSM_info::Ye21:
      Ye(2,1) = x;
      break;
   case CSE6SSM_info::Ye22:
      Ye(2,2) = x;
      break;
   case CSE6SSM_info::SigmaL:
      SigmaL = x;
      break;
   case CSE6SSM_info::KappaPr:
      KappaPr = x;
      break;
   case CSE6SSM_info::Sigmax:
      Sigmax = x;
      break;
   case CSE6SSM_info::gD00:
      gD(0,0) = x;
      break;
   case CSE6SSM_info::gD01:
      gD(0,1) = x;
      break;
   case CSE6SSM_info::gD02:
      gD(0,2) = x;
      break;
   case CSE6SSM_info::gD10:
      gD(1,0) = x;
      break;
   case CSE6SSM_info::gD11:
      gD(1,1) = x;
      break;
   case CSE6SSM_info::gD12:
      gD(1,2) = x;
      break;
   case CSE6SSM_info::gD20:
      gD(2,0) = x;
      break;
   case CSE6SSM_info::gD21:
      gD(2,1) = x;
      break;
   case CSE6SSM_info::gD22:
      gD(2,2) = x;
      break;
   case CSE6SSM_info::Kappa00:
      Kappa(0,0) = x;
      break;
   case CSE6SSM_info::Kappa01:
      Kappa(0,1) = x;
      break;
   case CSE6SSM_info::Kappa02:
      Kappa(0,2) = x;
      break;
   case CSE6SSM_info::Kappa10:
      Kappa(1,0) = x;
      break;
   case CSE6SSM_info::Kappa11:
      Kappa(1,1) = x;
      break;
   case CSE6SSM_info::Kappa12:
      Kappa(1,2) = x;
      break;
   case CSE6SSM_info::Kappa20:
      Kappa(2,0) = x;
      break;
   case CSE6SSM_info::Kappa21:
      Kappa(2,1) = x;
      break;
   case CSE6SSM_info::Kappa22:
      Kappa(2,2) = x;
      break;
   case CSE6SSM_info::Lambda1200:
      Lambda12(0,0) = x;
      break;
   case CSE6SSM_info::Lambda1201:
      Lambda12(0,1) = x;
      break;
   case CSE6SSM_info::Lambda1210:
      Lambda12(1,0) = x;
      break;
   case CSE6SSM_info::Lambda1211:
      Lambda12(1,1) = x;
      break;
   case CSE6SSM_info::Lambdax:
      Lambdax = x;
      break;
   case CSE6SSM_info::fu00:
      fu(0,0) = x;
      break;
   case CSE6SSM_info::fu01:
      fu(0,1) = x;
      break;
   case CSE6SSM_info::fu10:
      fu(1,0) = x;
      break;
   case CSE6SSM_info::fu11:
      fu(1,1) = x;
      break;
   case CSE6SSM_info::fu20:
      fu(2,0) = x;
      break;
   case CSE6SSM_info::fu21:
      fu(2,1) = x;
      break;
   case CSE6SSM_info::fd00:
      fd(0,0) = x;
      break;
   case CSE6SSM_info::fd01:
      fd(0,1) = x;
      break;
   case CSE6SSM_info::fd10:
      fd(1,0) = x;
      break;
   case CSE6SSM_info::fd11:
      fd(1,1) = x;
      break;
   case CSE6SSM_info::fd20:
      fd(2,0) = x;
      break;
   case CSE6SSM_info::fd21:
      fd(2,1) = x;
      break;
   case CSE6SSM_info::Yu00:
      Yu(0,0) = x;
      break;
   case CSE6SSM_info::Yu01:
      Yu(0,1) = x;
      break;
   case CSE6SSM_info::Yu02:
      Yu(0,2) = x;
      break;
   case CSE6SSM_info::Yu10:
      Yu(1,0) = x;
      break;
   case CSE6SSM_info::Yu11:
      Yu(1,1) = x;
      break;
   case CSE6SSM_info::Yu12:
      Yu(1,2) = x;
      break;
   case CSE6SSM_info::Yu20:
      Yu(2,0) = x;
      break;
   case CSE6SSM_info::Yu21:
      Yu(2,1) = x;
      break;
   case CSE6SSM_info::Yu22:
      Yu(2,2) = x;
      break;
   case CSE6SSM_info::MuPr:
      MuPr = x;
      break;
   case CSE6SSM_info::MuPhi:
      MuPhi = x;
      break;
   case CSE6SSM_info::XiF:
      XiF = x;
      break;
   case CSE6SSM_info::g1:
      g1 = x;
      break;
   case CSE6SSM_info::g2:
      g2 = x;
      break;
   case CSE6SSM_info::g3:
      g3 = x;
      break;
   case CSE6SSM_info::g1p:
      g1p = x;
      break;
   case CSE6SSM_info::vd:
      vd = x;
      break;
   case CSE6SSM_info::vu:
      vu = x;
      break;
   case CSE6SSM_info::vs:
      vs = x;
      break;
   case CSE6SSM_info::vsb:
      vsb = x;
      break;
   case CSE6SSM_info::vphi:
      vphi = x;
      break;
   case CSE6SSM_info::QS:
      QS = x;
      break;
   case CSE6SSM_info::TYd00:
      TYd(0,0) = x;
      break;
   case CSE6SSM_info::TYd01:
      TYd(0,1) = x;
      break;
   case CSE6SSM_info::TYd02:
      TYd(0,2) = x;
      break;
   case CSE6SSM_info::TYd10:
      TYd(1,0) = x;
      break;
   case CSE6SSM_info::TYd11:
      TYd(1,1) = x;
      break;
   case CSE6SSM_info::TYd12:
      TYd(1,2) = x;
      break; 
   case CSE6SSM_info::TYd20:
      TYd(2,0) = x;
      break;
   case CSE6SSM_info::TYd21:
      TYd(2,1) = x;
      break;
   case CSE6SSM_info::TYd22:
      TYd(2,2) = x;
      break;
   case CSE6SSM_info::ThE00:
      ThE(0,0) = x;
      break;
   case CSE6SSM_info::ThE01:
      ThE(0,1) = x;
      break;
   case CSE6SSM_info::ThE10:
      ThE(1,0) = x;
      break;
   case CSE6SSM_info::ThE11:
      ThE(1,1) = x;
      break;
   case CSE6SSM_info::ThE20:
      ThE(2,0) = x;
      break;
   case CSE6SSM_info::ThE21:
      ThE(2,1) = x;
      break;
   case CSE6SSM_info::TYe00:
      TYe(0,0) = x;
      break;
   case CSE6SSM_info::TYe01:
      TYe(0,1) = x;
      break;
   case CSE6SSM_info::TYe02:
      TYe(0,2) = x;
      break;
   case CSE6SSM_info::TYe10:
      TYe(1,0) = x;
      break;
   case CSE6SSM_info::TYe11:
      TYe(1,1) = x;
      break;
   case CSE6SSM_info::TYe12:
      TYe(1,2) = x;
      break;
   case CSE6SSM_info::TYe20:
      TYe(2,0) = x;
      break;
   case CSE6SSM_info::TYe21:
      TYe(2,1) = x;
      break;
   case CSE6SSM_info::TYe22:
      TYe(2,2) = x;
      break;
   case CSE6SSM_info::TSigmaL:
      TSigmaL = x;
      break;
   case CSE6SSM_info::TKappaPr:
      TKappaPr = x;
      break;
   case CSE6SSM_info::TSigmax:
      TSigmax = x;
      break;
   case CSE6SSM_info::TgD00:
      TgD(0,0) = x;
      break;
   case CSE6SSM_info::TgD01:
      TgD(0,1) = x;
      break;
   case CSE6SSM_info::TgD02:
      TgD(0,2) = x;
      break;
   case CSE6SSM_info::TgD10:
      TgD(1,0) = x;
      break;
   case CSE6SSM_info::TgD11:
      TgD(1,1) = x;
      break;
   case CSE6SSM_info::TgD12:
      TgD(1,2) = x;
      break;
   case CSE6SSM_info::TgD20:
      TgD(2,0) = x;
      break;
   case CSE6SSM_info::TgD21:
      TgD(2,1) = x;
      break;
   case CSE6SSM_info::TgD22:
      TgD(2,2) = x;
      break;
   case CSE6SSM_info::TKappa00:
      TKappa(0,0) = x;
      break;
   case CSE6SSM_info::TKappa01:
      TKappa(0,1) = x;
      break;
   case CSE6SSM_info::TKappa02:
      TKappa(0,2) = x;
      break;
   case CSE6SSM_info::TKappa10:
      TKappa(1,0) = x;
      break;
   case CSE6SSM_info::TKappa11:
      TKappa(1,1) = x;
      break;
   case CSE6SSM_info::TKappa12:
      TKappa(1,2) = x;
      break;
   case CSE6SSM_info::TKappa20:
      TKappa(2,0) = x;
      break;
   case CSE6SSM_info::TKappa21:
      TKappa(2,1) = x;
      break;
   case CSE6SSM_info::TKappa22:
      TKappa(2,2) = x;
      break;
   case CSE6SSM_info::TLambda1200:
      TLambda12(0,0) = x;
      break;
   case CSE6SSM_info::TLambda1201:
      TLambda12(0,1) = x;
      break;
   case CSE6SSM_info::TLambda1210:
      TLambda12(1,0) = x;
      break;
   case CSE6SSM_info::TLambda1211:
      TLambda12(1,1) = x;
      break;
   case CSE6SSM_info::TLambdax:
      TLambdax = x;
      break;
   case CSE6SSM_info::Tfu00:
      Tfu(0,0) = x;
      break;
   case CSE6SSM_info::Tfu01:
      Tfu(0,1) = x;
      break;
   case CSE6SSM_info::Tfu10:
      Tfu(1,0) = x;
      break;
   case CSE6SSM_info::Tfu11:
      Tfu(1,1) = x;
      break;
   case CSE6SSM_info::Tfu20:
      Tfu(2,0) = x;
      break;
   case CSE6SSM_info::Tfu21:
      Tfu(2,1) = x;
      break;
   case CSE6SSM_info::Tfd00:
      Tfd(0,0) = x;
      break;
   case CSE6SSM_info::Tfd01:
      Tfd(0,1) = x;
      break;
   case CSE6SSM_info::Tfd10:
      Tfd(1,0) = x;
      break;
   case CSE6SSM_info::Tfd11:
      Tfd(1,1) = x;
      break;
   case CSE6SSM_info::Tfd20:
      Tfd(2,0) = x;
      break;
   case CSE6SSM_info::Tfd21:
      Tfd(2,1) = x;
      break;
   case CSE6SSM_info::TYu00:
      TYu(0,0) = x;
      break;
   case CSE6SSM_info::TYu01:
      TYu(0,1) = x;
      break;
   case CSE6SSM_info::TYu02:
      TYu(0,2) = x;
      break;
   case CSE6SSM_info::TYu10:
      TYu(1,0) = x;
      break;
   case CSE6SSM_info::TYu11:
      TYu(1,1) = x;
      break;
   case CSE6SSM_info::TYu12:
      TYu(1,2) = x;
      break;
   case CSE6SSM_info::TYu20:
      TYu(2,0) = x;
      break;
   case CSE6SSM_info::TYu21:
      TYu(2,1) = x;
      break;
   case CSE6SSM_info::TYu22:
      TYu(2,2) = x;
      break;
   case CSE6SSM_info::BMuPr:
      BMuPr = x;
      break;
   case CSE6SSM_info::BMuPhi:
      BMuPhi = x;
      break;
   case CSE6SSM_info::LXiF:
      LXiF = x;
      break;
   case CSE6SSM_info::mq200:
      mq2(0,0) = x;
      break;
   case CSE6SSM_info::mq201:
      mq2(0,1) = x;
      break;
   case CSE6SSM_info::mq202:
      mq2(0,2) = x;
      break;
   case CSE6SSM_info::mq210:
      mq2(1,0) = x;
      break;
   case CSE6SSM_info::mq211:
      mq2(1,1) = x;
      break;
   case CSE6SSM_info::mq212:
      mq2(1,2) = x;
      break;
   case CSE6SSM_info::mq220:
      mq2(2,0) = x;
      break;
   case CSE6SSM_info::mq221:
      mq2(2,1) = x;
      break;
   case CSE6SSM_info::mq222:
      mq2(2,2) = x;
      break;
   case CSE6SSM_info::ml200:
      ml2(0,0) = x;
      break;
   case CSE6SSM_info::ml201:
      ml2(0,1) = x;
      break;
   case CSE6SSM_info::ml202:
      ml2(0,2) = x;
      break;
   case CSE6SSM_info::ml210:
      ml2(1,0) = x;
      break;
   case CSE6SSM_info::ml211:
      ml2(1,1) = x;
      break;
   case CSE6SSM_info::ml212:
      ml2(1,2) = x;
      break;
   case CSE6SSM_info::ml220:
      ml2(2,0) = x;
      break;
   case CSE6SSM_info::ml221:
      ml2(2,1) = x;
      break;
   case CSE6SSM_info::ml222:
      ml2(2,2) = x;
      break;
   case CSE6SSM_info::mHd2:
      mHd2 = x;
      break;
   case CSE6SSM_info::mHu2:
      mHu2 = x;
      break;
   case CSE6SSM_info::md200:
      md2(0,0) = x;
      break;
   case CSE6SSM_info::md201:
      md2(0,1) = x;
      break;
   case CSE6SSM_info::md202:
      md2(0,2) = x;
      break;
   case CSE6SSM_info::md210:
      md2(1,0) = x;
      break;
   case CSE6SSM_info::md211:
      md2(1,1) = x;
      break;
   case CSE6SSM_info::md212:
      md2(1,2) = x;
      break;
   case CSE6SSM_info::md220:
      md2(2,0) = x;
      break;
   case CSE6SSM_info::md221:
      md2(2,1) = x;
      break;
   case CSE6SSM_info::md222:
      md2(2,2) = x;
      break;
   case CSE6SSM_info::mu200:
      mu2(0,0) = x;
      break;
   case CSE6SSM_info::mu201:
      mu2(0,1) = x;
      break;
   case CSE6SSM_info::mu202:
      mu2(0,2) = x;
      break;
   case CSE6SSM_info::mu210:
      mu2(1,0) = x;
      break;
   case CSE6SSM_info::mu211:
      mu2(1,1) = x;
      break;
   case CSE6SSM_info::mu212:
      mu2(1,2) = x;
      break;
   case CSE6SSM_info::mu220:
      mu2(2,0) = x;
      break;
   case CSE6SSM_info::mu221:
      mu2(2,1) = x;
      break;
   case CSE6SSM_info::mu222:
      mu2(2,2) = x;
      break;
   case CSE6SSM_info::me200:
      me2(0,0) = x;
      break;
   case CSE6SSM_info::me201:
      me2(0,1) = x;
      break;
   case CSE6SSM_info::me202:
      me2(0,2) = x;
      break;
   case CSE6SSM_info::me210:
      me2(1,0) = x;
      break;
   case CSE6SSM_info::me211:
      me2(1,1) = x;
      break;
   case CSE6SSM_info::me212:
      me2(1,2) = x;
      break;
   case CSE6SSM_info::me220:
      me2(2,0) = x;
      break;
   case CSE6SSM_info::me221:
      me2(2,1) = x;
      break;
   case CSE6SSM_info::me222:
      me2(2,2) = x;
      break;
   case CSE6SSM_info::ms2:
      ms2 = x;
      break;
   case CSE6SSM_info::msbar2:
      msbar2 = x;
      break;
   case CSE6SSM_info::mH1I200:
      mH1I2(0,0) = x;
      break;
   case CSE6SSM_info::mH1I201:
      mH1I2(0,1) = x;
      break;
   case CSE6SSM_info::mH1I210:
      mH1I2(1,0) = x;
      break;
   case CSE6SSM_info::mH1I211:
      mH1I2(1,1) = x;
      break;
   case CSE6SSM_info::mH2I200:
      mH2I2(0,0) = x;
      break;
   case CSE6SSM_info::mH2I201:
      mH2I2(0,1) = x;
      break;
   case CSE6SSM_info::mH2I210:
      mH2I2(1,0) = x;
      break;
   case CSE6SSM_info::mH2I211:
      mH2I2(1,1) = x;
      break;
   case CSE6SSM_info::mSI200:
      mSI2(0,0) = x;
      break;
   case CSE6SSM_info::mSI201:
      mSI2(0,1) = x;
      break;
   case CSE6SSM_info::mSI202:
      mSI2(0,2) = x;
      break;
   case CSE6SSM_info::mSI210:
      mSI2(1,0) = x;
      break;
   case CSE6SSM_info::mSI211:
      mSI2(1,1) = x;
      break;
   case CSE6SSM_info::mSI212:
      mSI2(1,2) = x;
      break;
   case CSE6SSM_info::mSI220:
      mSI2(2,0) = x;
      break;
   case CSE6SSM_info::mSI221:
      mSI2(2,1) = x;
      break;
   case CSE6SSM_info::mSI222:
      mSI2(2,2) = x;
      break;
   case CSE6SSM_info::mDx200:
      mDx2(0,0) = x;
      break;
   case CSE6SSM_info::mDx201:
      mDx2(0,1) = x;
      break;
   case CSE6SSM_info::mDx202:
      mDx2(0,2) = x;
      break;
   case CSE6SSM_info::mDx210:
      mDx2(1,0) = x;
      break;
   case CSE6SSM_info::mDx211:
      mDx2(1,1) = x;
      break;
   case CSE6SSM_info::mDx212:
      mDx2(1,2) = x;
      break;
   case CSE6SSM_info::mDx220:
      mDx2(2,0) = x;
      break;
   case CSE6SSM_info::mDx221:
      mDx2(2,1) = x;
      break;
   case CSE6SSM_info::mDx222:
      mDx2(2,2) = x;
      break;
   case CSE6SSM_info::mDxbar200:
      mDxbar2(0,0) = x;
      break;
   case CSE6SSM_info::mDxbar201:
      mDxbar2(0,1) = x;
      break;
   case CSE6SSM_info::mDxbar202:
      mDxbar2(0,2) = x;
      break;
   case CSE6SSM_info::mDxbar210:
      mDxbar2(1,0) = x;
      break;
   case CSE6SSM_info::mDxbar211:
      mDxbar2(1,1) = x;
      break;
   case CSE6SSM_info::mDxbar212:
      mDxbar2(1,2) = x;
      break;
   case CSE6SSM_info::mDxbar220:
      mDxbar2(2,0) = x;
      break;
   case CSE6SSM_info::mDxbar221:
      mDxbar2(2,1) = x;
      break;
   case CSE6SSM_info::mDxbar222:
      mDxbar2(2,2) = x;
      break;
   case CSE6SSM_info::mHp2:
      mHp2 = x;
      break;
   case CSE6SSM_info::mHpbar2:
      mHpbar2 = x;
      break;
   case CSE6SSM_info::mphi2:
      mphi2 = x;
      break;
   case CSE6SSM_info::MassB:
      MassB = x;
      break;
   case CSE6SSM_info::MassWB:
      MassWB = x;
      break;
   case CSE6SSM_info::MassG:
      MassG = x;
      break;
   case CSE6SSM_info::MassBp:
      MassBp = x;
      break;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

std::ostream& operator<<(std::ostream& ostr, const CSE6SSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
