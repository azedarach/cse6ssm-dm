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
 * Calculates the one-loop beta function of SigmaL.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_SigmaL_one_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   double beta_SigmaL;

   beta_SigmaL = Re(0.2*oneOver16PiSqr*SigmaL*(15*tracegDAdjgD + 5*
      tracehEAdjhE + 10*AbsSqr(KappaPr) + 5*AbsSqr(Sigmax) + 20*AbsSqr(SigmaL)
      - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2)));


   return beta_SigmaL;
}

/**
 * Calculates the two-loop beta function of SigmaL.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_SigmaL_two_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracegDAdjgDgDAdjgD = TRACE_STRUCT.tracegDAdjgDgDAdjgD;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjhEhEAdjhE = TRACE_STRUCT.tracehEAdjhEhEAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;


   double beta_SigmaL;

   beta_SigmaL = Re(-0.02*twoLoop*SigmaL*(-297*Power(g1,4) - 192*Power(
      g1p,4) - 825*Power(g2,4) + 50*tracefuAdjhEhEAdjfu + 450*
      tracegDAdjgDgDAdjgD + 150*tracegDAdjgDTpYdconjYd + 150*
      tracegDAdjgDTpYuconjYu + 150*tracegDAdjKappaKappaAdjgD + 150*
      tracehEAdjhEhEAdjhE + 100*tracehEAdjhEYeAdjYe + 50*
      tracehEAdjLambda12Lambda12AdjhE + 450*tracegDAdjgD*AbsSqr(SigmaL) + 150*
      tracehEAdjhE*AbsSqr(SigmaL) + 200*AbsSqr(KappaPr)*(AbsSqr(Sigmax) + 3*
      AbsSqr(SigmaL)) + 20*tracegDAdjgD*Sqr(g1) - 60*tracehEAdjhE*Sqr(g1) - 60*
      AbsSqr(SigmaL)*Sqr(g1) - 45*tracegDAdjgD*Sqr(g1p) - 15*tracehEAdjhE*Sqr(
      g1p) - 40*AbsSqr(SigmaL)*Sqr(g1p) - 36*Sqr(g1)*Sqr(g1p) - 300*AbsSqr(
      SigmaL)*Sqr(g2) - 90*Sqr(g1)*Sqr(g2) - 60*Sqr(g1p)*Sqr(g2) - 800*
      tracegDAdjgD*Sqr(g3) - Power(g1p,4)*Sqr(QS) - 5*AbsSqr(Sigmax)*(-30*
      traceKappaAdjKappa - 20*traceLambda12AdjLambda12 - 20*AbsSqr(Lambdax) -
      20*AbsSqr(SigmaL) + Sqr(g1p)*Sqr(QS)) + 400*Sqr(Conj(KappaPr))*Sqr(
      KappaPr) + 100*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 500*Sqr(Conj(SigmaL))*Sqr(
      SigmaL)));


   return beta_SigmaL;
}

} // namespace flexiblesusy
