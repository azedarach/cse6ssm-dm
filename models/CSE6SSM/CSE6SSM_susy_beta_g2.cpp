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
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(4*Power(g2,3)*oneOver16PiSqr);


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_g2;

   beta_g2 = Re(0.2*Power(g2,3)*twoLoop*(-10*tracefdAdjfd - 10*
      tracefuAdjfu - 30*tracegDAdjgD - 10*tracehEAdjhE - 10*
      traceLambda12AdjLambda12 - 30*traceYdAdjYd - 10*traceYeAdjYe - 30*
      traceYuAdjYu - 10*AbsSqr(Lambdax) - 10*AbsSqr(SigmaL) + 18*Sqr(g1) + 17*
      Sqr(g1p) + 230*Sqr(g2) + 120*Sqr(g3)));


   return beta_g2;
}

} // namespace flexiblesusy
