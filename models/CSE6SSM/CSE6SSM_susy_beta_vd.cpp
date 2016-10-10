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
 * Calculates the one-loop beta function of vd.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vd_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = Re(0.05*oneOver16PiSqr*vd*(-20*tracefdAdjfd - 60*
      traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(Lambdax) + 6*Sqr(g1) + 9*Sqr(
      g1p) + 30*Sqr(g2)));


   return beta_vd;
}

/**
 * Calculates the two-loop beta function of vd.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vd_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfdfdAdjfd = TRACE_STRUCT.tracefdAdjfdfdAdjfd;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;


   double beta_vd;

   beta_vd = Re(-0.00125*twoLoop*vd*(1188*Power(g1,4) + 1773*Power(g1p,4)
      + 2900*Power(g2,4) - 2400*tracefdAdjfdfdAdjfd - 1600*tracefdAdjfdfuAdjfu
      - 2400*tracegDAdjgDTpYdconjYd - 1600*tracehEAdjhEYeAdjYe - 800*
      traceLambda12AdjLambda12Tpfdconjfd - 7200*traceYdAdjYdYdAdjYd - 2400*
      traceYdAdjYuYuAdjYd - 2400*traceYeAdjYeYeAdjYe + 400*traceYdAdjYd*Sqr(g1)
      + 1200*traceYeAdjYe*Sqr(g1) + 600*traceYdAdjYd*Sqr(g1p) + 200*
      traceYeAdjYe*Sqr(g1p) - 36*Sqr(g1)*Sqr(g1p) + 3600*traceYdAdjYd*Sqr(g2) +
      1200*traceYeAdjYe*Sqr(g2) + 360*Sqr(g1)*Sqr(g2) + 540*Sqr(g1p)*Sqr(g2) +
      40*tracefdAdjfd*(6*Sqr(g1) + 29*Sqr(g1p) + 30*Sqr(g2)) + 12800*
      traceYdAdjYd*Sqr(g3) + 9*Power(g1p,4)*Sqr(QS) + 40*AbsSqr(Lambdax)*(-20*
      tracefuAdjfu - 60*traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 60*
      traceYuAdjYu - 20*AbsSqr(Sigmax) + 6*Sqr(g1) + 4*Sqr(g1p) + 30*Sqr(g2) +
      Sqr(g1p)*Sqr(QS)) - 2400*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vd;
}

} // namespace flexiblesusy
