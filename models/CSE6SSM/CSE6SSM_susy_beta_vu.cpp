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

// File generated at Wed 3 Jun 2015 23:42:49

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(-10*tracefuAdjfu - 30*traceYuAdjYu
      - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 2*Sqr(g1p) + 15*Sqr(g2)));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(297*Power(g1,4) + 192*Power(g1p,4) +
      725*Power(g2,4) - 400*tracefdAdjfdfuAdjfu - 600*tracefuAdjfufuAdjfu - 200
      *tracefuAdjhEhEAdjfu - 200*tracefuAdjLambda12Lambda12Adjfu - 600*
      tracegDAdjgDTpYuconjYu - 600*traceYdAdjYuYuAdjYd - 1800*
      traceYuAdjYuYuAdjYu + 340*traceYuAdjYu*Sqr(g1) + 60*traceYuAdjYu*Sqr(g1p)
      + 36*Sqr(g1)*Sqr(g1p) + 900*traceYuAdjYu*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) +
      60*Sqr(g1p)*Sqr(g2) + 20*tracefuAdjfu*(3*Sqr(g1) + 17*Sqr(g1p) + 15*Sqr(
      g2)) + 3200*traceYuAdjYu*Sqr(g3) + Power(g1p,4)*Sqr(QS) + 10*AbsSqr(
      Lambdax)*(-20*tracefdAdjfd - 60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 60*traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(
      Sigmax) + 6*Sqr(g1) + 9*Sqr(g1p) + 30*Sqr(g2) + Sqr(g1p)*Sqr(QS)) - 600*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vu;
}

} // namespace flexiblesusy
