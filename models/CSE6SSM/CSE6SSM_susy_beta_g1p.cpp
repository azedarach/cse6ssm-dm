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
 * Calculates the one-loop beta function of g1p.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g1p_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g1p;

   beta_g1p = Re(0.05*Power(g1p,3)*oneOver16PiSqr*(188 + Sqr(QS)));


   return beta_g1p;
}

/**
 * Calculates the two-loop beta function of g1p.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_g1p_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_g1p;

   beta_g1p = Re(0.005*Power(g1p,3)*twoLoop*(-760*tracefdAdjfd - 760*
      tracefuAdjfu - 840*tracegDAdjgD - 280*tracehEAdjhE - 390*
      traceKappaAdjKappa - 260*traceLambda12AdjLambda12 - 840*traceYdAdjYd -
      280*traceYeAdjYe - 360*traceYuAdjYu - 160*AbsSqr(SigmaL) + 648*Sqr(g1) +
      1832*Sqr(g1p) + Power(QS,4)*Sqr(g1p) + 2040*Sqr(g2) + 4800*Sqr(g3) - 30*
      traceKappaAdjKappa*Sqr(QS) - 20*traceLambda12AdjLambda12*Sqr(QS) - 20*
      AbsSqr(Sigmax)*Sqr(QS) - 20*AbsSqr(Lambdax)*(13 + Sqr(QS))));


   return beta_g1p;
}

} // namespace flexiblesusy
