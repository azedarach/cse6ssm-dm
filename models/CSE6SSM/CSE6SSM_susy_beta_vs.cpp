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
 * Calculates the one-loop beta function of vs.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vs_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_vs;

   beta_vs = Re(oneOver16PiSqr*(-3*traceKappaAdjKappa*vs - 2*
      traceLambda12AdjLambda12*vs - 2*vs*AbsSqr(Lambdax) - vs*AbsSqr(Sigmax) +
      0.05*vs*Sqr(g1p)*Sqr(QS)));


   return beta_vs;
}

/**
 * Calculates the two-loop beta function of vs.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vs_two_loop(const Susy_traces& susy_traces) const
{
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


   double beta_vs;

   beta_vs = Re(-0.0025*twoLoop*vs*(Power(g1p,4)*Power(QS,4) - 800*
      tracefuAdjLambda12Lambda12Adjfu - 2400*tracegDAdjKappaKappaAdjgD - 800*
      tracehEAdjLambda12Lambda12AdjhE - 2400*traceKappaAdjKappaKappaAdjKappa -
      1600*traceLambda12AdjLambda12Lambda12AdjLambda12 - 800*
      traceLambda12AdjLambda12Tpfdconjfd + 320*traceKappaAdjKappa*Sqr(g1) + 480
      *traceLambda12AdjLambda12*Sqr(g1) + 780*traceKappaAdjKappa*Sqr(g1p) + 520
      *traceLambda12AdjLambda12*Sqr(g1p) + 2400*traceLambda12AdjLambda12*Sqr(g2
      ) + 40*AbsSqr(Lambdax)*(-20*tracefdAdjfd - 20*tracefuAdjfu - 60*
      traceYdAdjYd - 20*traceYeAdjYe - 60*traceYuAdjYu + 12*Sqr(g1) + 13*Sqr(
      g1p) + 60*Sqr(g2)) + 6400*traceKappaAdjKappa*Sqr(g3) + 94*Power(g1p,4)*
      Sqr(QS) + 20*AbsSqr(Sigmax)*(-40*AbsSqr(KappaPr) - 40*AbsSqr(SigmaL) +
      Sqr(g1p)*Sqr(QS)) - 1600*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 800*Sqr(Conj(
      Sigmax))*Sqr(Sigmax)));


   return beta_vs;
}

} // namespace flexiblesusy
