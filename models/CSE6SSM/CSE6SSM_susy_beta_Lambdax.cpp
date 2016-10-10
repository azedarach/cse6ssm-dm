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

// File generated at Wed 3 Jun 2015 23:42:46

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = Re(0.05*oneOver16PiSqr*Lambdax*(20*tracefdAdjfd + 20*
      tracefuAdjfu + 60*traceKappaAdjKappa + 40*traceLambda12AdjLambda12 + 60*
      traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu + 80*AbsSqr(Lambdax) +
      20*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) - Sqr(g1p)*Sqr(
      QS)));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfdfdAdjfd = TRACE_STRUCT.tracefdAdjfdfdAdjfd;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.0025*twoLoop*Lambdax*(-2376*Power(g1,4) - 2541*
      Power(g1p,4) - 6600*Power(g2,4) - 2*Power(g1p,4)*Power(QS,4) + 1200*
      tracefdAdjfdfdAdjfd + 1600*tracefdAdjfdfuAdjfu + 1200*tracefuAdjfufuAdjfu
      + 400*tracefuAdjhEhEAdjfu + 1200*tracefuAdjLambda12Lambda12Adjfu + 1200*
      tracegDAdjgDTpYdconjYd + 1200*tracegDAdjgDTpYuconjYu + 2400*
      tracegDAdjKappaKappaAdjgD + 800*tracehEAdjhEYeAdjYe + 800*
      tracehEAdjLambda12Lambda12AdjhE + 2400*traceKappaAdjKappaKappaAdjKappa +
      1600*traceLambda12AdjLambda12Lambda12AdjLambda12 + 1200*
      traceLambda12AdjLambda12Tpfdconjfd + 3600*traceYdAdjYdYdAdjYd + 2400*
      traceYdAdjYuYuAdjYd + 1200*traceYeAdjYeYeAdjYe + 3600*traceYuAdjYuYuAdjYu
      + 800*AbsSqr(KappaPr)*AbsSqr(Sigmax) + 800*AbsSqr(Sigmax)*AbsSqr(SigmaL)
      - 320*traceKappaAdjKappa*Sqr(g1) - 480*traceLambda12AdjLambda12*Sqr(g1)
      + 160*traceYdAdjYd*Sqr(g1) - 480*traceYeAdjYe*Sqr(g1) - 320*traceYuAdjYu*
      Sqr(g1) - 400*tracefdAdjfd*Sqr(g1p) - 600*tracefuAdjfu*Sqr(g1p) - 780*
      traceKappaAdjKappa*Sqr(g1p) - 520*traceLambda12AdjLambda12*Sqr(g1p) + 240
      *traceYdAdjYd*Sqr(g1p) + 80*traceYeAdjYe*Sqr(g1p) + 120*traceYuAdjYu*Sqr(
      g1p) - 108*Sqr(g1)*Sqr(g1p) - 2400*traceLambda12AdjLambda12*Sqr(g2) - 720
      *Sqr(g1)*Sqr(g2) - 780*Sqr(g1p)*Sqr(g2) - 40*AbsSqr(Lambdax)*(-30*
      tracefdAdjfd - 30*tracefuAdjfu - 60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 90*traceYdAdjYd - 30*traceYeAdjYe - 90*
      traceYuAdjYu - 20*AbsSqr(Sigmax) + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2))
      - 6400*traceKappaAdjKappa*Sqr(g3) - 6400*traceYdAdjYd*Sqr(g3) - 6400*
      traceYuAdjYu*Sqr(g3) - 201*Power(g1p,4)*Sqr(QS) + 60*traceKappaAdjKappa*
      Sqr(g1p)*Sqr(QS) + 40*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS) + 4000*
      Sqr(Conj(Lambdax))*Sqr(Lambdax) + 800*Sqr(Conj(Sigmax))*Sqr(Sigmax)));


   return beta_Lambdax;
}

} // namespace flexiblesusy
