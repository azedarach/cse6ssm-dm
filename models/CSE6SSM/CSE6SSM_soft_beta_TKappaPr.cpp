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

// File generated at Wed 3 Jun 2015 23:43:02

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappaPr.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TKappaPr_one_loop(const Soft_traces& soft_traces) const
{


   double beta_TKappaPr;

   beta_TKappaPr = Re(3*oneOver16PiSqr*(6*AbsSqr(KappaPr)*TKappaPr + Conj
      (Sigmax)*(Sigmax*TKappaPr + 2*KappaPr*TSigmax) + 2*Conj(SigmaL)*(SigmaL*
      TKappaPr + 2*KappaPr*TSigmaL)));


   return beta_TKappaPr;
}

/**
 * Calculates the two-loop beta function of TKappaPr.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TKappaPr_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_TKappaPr;

   beta_TKappaPr = Re(-0.3*twoLoop*(20*Sigmax*Sqr(Conj(Sigmax))*(Sigmax*
      TKappaPr + 4*KappaPr*TSigmax) + Conj(Sigmax)*(Sigmax*(30*
      traceKappaAdjKappa + 20*traceLambda12AdjLambda12 + 120*AbsSqr(KappaPr) +
      20*AbsSqr(Lambdax) - Sqr(g1p)*Sqr(QS))*TKappaPr + 2*KappaPr*(Sigmax*(30*
      traceAdjKappaTKappa + 20*traceAdjLambda12TLambda12 + MassBp*Sqr(g1p)*Sqr(
      QS)) + (30*traceKappaAdjKappa + 20*traceLambda12AdjLambda12 + 40*AbsSqr(
      KappaPr) - Sqr(g1p)*Sqr(QS))*TSigmax + 20*Conj(Lambdax)*(Sigmax*TLambdax
      + Lambdax*TSigmax))) + 4*(100*Sqr(Conj(KappaPr))*Sqr(KappaPr)*TKappaPr +
      10*SigmaL*Sqr(Conj(SigmaL))*(SigmaL*TKappaPr + 4*KappaPr*TSigmaL) + Conj(
      SigmaL)*(SigmaL*(15*tracegDAdjgD + 5*tracehEAdjhE + 60*AbsSqr(KappaPr) -
      3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2))*TKappaPr + 2*KappaPr*(SigmaL*(15*
      traceAdjgDTgD + 5*traceAdjhEThE + 3*MassB*Sqr(g1) + 2*MassBp*Sqr(g1p) +
      15*MassWB*Sqr(g2)) + (15*tracegDAdjgD + 5*tracehEAdjhE + 20*AbsSqr(
      KappaPr) - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2))*TSigmaL)))));


   return beta_TKappaPr;
}

} // namespace flexiblesusy
