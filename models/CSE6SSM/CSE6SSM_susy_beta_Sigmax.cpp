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

// File generated at Wed 3 Jun 2015 23:42:44

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Sigmax.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_Sigmax_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_Sigmax;

   beta_Sigmax = Re(oneOver16PiSqr*(3*traceKappaAdjKappa*Sigmax + 2*
      traceLambda12AdjLambda12*Sigmax + 2*AbsSqr(KappaPr)*Sigmax + 2*AbsSqr(
      Lambdax)*Sigmax + 2*AbsSqr(SigmaL)*Sigmax - 0.1*Sigmax*Sqr(g1p)*Sqr(QS) +
      3*Conj(Sigmax)*Sqr(Sigmax)));


   return beta_Sigmax;
}

/**
 * Calculates the two-loop beta function of Sigmax.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_Sigmax_two_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;


   double beta_Sigmax;

   beta_Sigmax = Re(-0.01*twoLoop*Sigmax*(-(Power(g1p,4)*Power(QS,4)) +
      200*tracefuAdjLambda12Lambda12Adjfu + 600*tracegDAdjKappaKappaAdjgD + 200
      *tracehEAdjLambda12Lambda12AdjhE + 600*traceKappaAdjKappaKappaAdjKappa +
      400*traceLambda12AdjLambda12Lambda12AdjLambda12 + 200*
      traceLambda12AdjLambda12Tpfdconjfd + 600*traceKappaAdjKappa*AbsSqr(Sigmax
      ) + 400*traceLambda12AdjLambda12*AbsSqr(Sigmax) + 600*tracegDAdjgD*AbsSqr
      (SigmaL) + 200*tracehEAdjhE*AbsSqr(SigmaL) + 400*AbsSqr(Sigmax)*AbsSqr(
      SigmaL) + 800*AbsSqr(KappaPr)*(AbsSqr(Sigmax) + AbsSqr(SigmaL)) - 80*
      traceKappaAdjKappa*Sqr(g1) - 120*traceLambda12AdjLambda12*Sqr(g1) - 120*
      AbsSqr(SigmaL)*Sqr(g1) - 195*traceKappaAdjKappa*Sqr(g1p) - 130*
      traceLambda12AdjLambda12*Sqr(g1p) - 80*AbsSqr(SigmaL)*Sqr(g1p) - 600*
      traceLambda12AdjLambda12*Sqr(g2) - 600*AbsSqr(SigmaL)*Sqr(g2) - 1600*
      traceKappaAdjKappa*Sqr(g3) - 94*Power(g1p,4)*Sqr(QS) + 15*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) + 10*traceLambda12AdjLambda12*Sqr(g1p
      )*Sqr(QS) - 10*AbsSqr(Sigmax)*Sqr(g1p)*Sqr(QS) - 10*AbsSqr(Lambdax)*(-20*
      tracefdAdjfd - 20*tracefuAdjfu - 60*traceYdAdjYd - 20*traceYeAdjYe - 60*
      traceYuAdjYu - 40*AbsSqr(Sigmax) + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2)
      - Sqr(g1p)*Sqr(QS)) + 800*Sqr(Conj(KappaPr))*Sqr(KappaPr) + 400*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 600*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 400*Sqr(Conj
      (SigmaL))*Sqr(SigmaL)));


   return beta_Sigmax;
}

} // namespace flexiblesusy
