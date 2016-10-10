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

// File generated at Wed 3 Jun 2015 23:42:48

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g3.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g3_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g3;

   beta_g3 = Re(0);


   return beta_g3;
}

/**
 * Calculates the two-loop beta function of g3.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g3_two_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_g3;

   beta_g3 = Re(Power(g3,3)*twoLoop*(-4*tracegDAdjgD - 2*
      traceKappaAdjKappa - 4*traceYdAdjYd - 4*traceYuAdjYu + 3*Sqr(g1) + 3*Sqr(
      g1p) + 9*Sqr(g2) + 48*Sqr(g3)));


   return beta_g3;
}

} // namespace flexiblesusy
