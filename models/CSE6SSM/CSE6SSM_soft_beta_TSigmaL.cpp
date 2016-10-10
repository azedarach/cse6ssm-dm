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
 * Calculates the one-loop beta function of TSigmaL.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TSigmaL_one_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;


   double beta_TSigmaL;

   beta_TSigmaL = Re(oneOver16PiSqr*(6*traceAdjgDTgD*SigmaL + 2*
      traceAdjhEThE*SigmaL + 1.2*MassB*SigmaL*Sqr(g1) + 0.8*MassBp*SigmaL*Sqr(
      g1p) + 6*MassWB*SigmaL*Sqr(g2) + 3*tracegDAdjgD*TSigmaL + tracehEAdjhE*
      TSigmaL + 12*AbsSqr(SigmaL)*TSigmaL - 0.6*Sqr(g1)*TSigmaL - 0.4*Sqr(g1p)*
      TSigmaL - 3*Sqr(g2)*TSigmaL + 2*Conj(KappaPr)*(2*SigmaL*TKappaPr +
      KappaPr*TSigmaL) + Conj(Sigmax)*(2*SigmaL*TSigmax + Sigmax*TSigmaL)));


   return beta_TSigmaL;
}

/**
 * Calculates the two-loop beta function of TSigmaL.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TSigmaL_two_loop(const Soft_traces& soft_traces) const
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
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjhEThEAdjfu = TRACE_STRUCT.tracefuAdjhEThEAdjfu;
   const double tracegDAdjgDgDAdjgD = TRACE_STRUCT.tracegDAdjgDgDAdjgD;
   const double tracegDAdjgDTgDAdjgD = TRACE_STRUCT.tracegDAdjgDTgDAdjgD;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjfuTfuAdjhE = TRACE_STRUCT.tracehEAdjfuTfuAdjhE;
   const double tracehEAdjhEhEAdjhE = TRACE_STRUCT.tracehEAdjhEhEAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjhEThEAdjhE = TRACE_STRUCT.tracehEAdjhEThEAdjhE;
   const double tracehEAdjhETYeAdjYe = TRACE_STRUCT.tracehEAdjhETYeAdjYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE;
   const double traceYeAdjYeThEAdjhE = TRACE_STRUCT.traceYeAdjYeThEAdjhE;
   const double traceKappaAdjgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjKappa;
   const double traceLambda12AdjhEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjLambda12;
   const double traceAdjgDTpYdconjYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjYdTgD;
   const double traceAdjgDTpYuconjYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjYuTgD;
   const double traceAdjYdTYdconjgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjgDTpgD;
   const double traceAdjYuTYuconjgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjgDTpgD;


   double beta_TSigmaL;

   beta_TSigmaL = Re(twoLoop*(-23.76*Power(g1,4)*MassB*SigmaL - 15.36*
      Power(g1p,4)*MassBp*SigmaL - 66*Power(g2,4)*MassWB*SigmaL - 6*
      traceAdjgDTpYdconjYdTgD*SigmaL - 6*traceAdjgDTpYuconjYuTgD*SigmaL - 6*
      traceAdjYdTYdconjgDTpgD*SigmaL - 6*traceAdjYuTYuconjgDTpgD*SigmaL - 2*
      tracefuAdjhEThEAdjfu*SigmaL - 36*tracegDAdjgDTgDAdjgD*SigmaL - 6*
      tracegDAdjKappaTKappaAdjgD*SigmaL - 2*tracehEAdjfuTfuAdjhE*SigmaL - 12*
      tracehEAdjhEThEAdjhE*SigmaL - 4*tracehEAdjhETYeAdjYe*SigmaL - 2*
      tracehEAdjLambda12TLambda12AdjhE*SigmaL - 6*traceKappaAdjgDTgDAdjKappa*
      SigmaL - 2*traceLambda12AdjhEThEAdjLambda12*SigmaL - 4*
      traceYeAdjYeThEAdjhE*SigmaL - 0.8*traceAdjgDTgD*SigmaL*Sqr(g1) + 2.4*
      traceAdjhEThE*SigmaL*Sqr(g1) + 0.8*MassB*tracegDAdjgD*SigmaL*Sqr(g1) -
      2.4*MassB*tracehEAdjhE*SigmaL*Sqr(g1) + 1.8*traceAdjgDTgD*SigmaL*Sqr(g1p)
      + 0.6*traceAdjhEThE*SigmaL*Sqr(g1p) - 1.8*MassBp*tracegDAdjgD*SigmaL*Sqr
      (g1p) - 0.6*MassBp*tracehEAdjhE*SigmaL*Sqr(g1p) - 1.44*MassB*SigmaL*Sqr(
      g1)*Sqr(g1p) - 1.44*MassBp*SigmaL*Sqr(g1)*Sqr(g1p) - 3.6*MassB*SigmaL*Sqr
      (g1)*Sqr(g2) - 3.6*MassWB*SigmaL*Sqr(g1)*Sqr(g2) - 2.4*MassBp*SigmaL*Sqr(
      g1p)*Sqr(g2) - 2.4*MassWB*SigmaL*Sqr(g1p)*Sqr(g2) + 32*traceAdjgDTgD*
      SigmaL*Sqr(g3) - 32*MassG*tracegDAdjgD*SigmaL*Sqr(g3) - 0.08*Power(g1p,4)
      *MassBp*SigmaL*Sqr(QS) - 32*KappaPr*SigmaL*Sqr(Conj(KappaPr))*TKappaPr +
      5.94*Power(g1,4)*TSigmaL + 3.84*Power(g1p,4)*TSigmaL + 16.5*Power(g2,4)*
      TSigmaL - tracefuAdjhEhEAdjfu*TSigmaL - 9*tracegDAdjgDgDAdjgD*TSigmaL - 3
      *tracegDAdjgDTpYdconjYd*TSigmaL - 3*tracegDAdjgDTpYuconjYu*TSigmaL - 3*
      tracegDAdjKappaKappaAdjgD*TSigmaL - 3*tracehEAdjhEhEAdjhE*TSigmaL - 2*
      tracehEAdjhEYeAdjYe*TSigmaL - tracehEAdjLambda12Lambda12AdjhE*TSigmaL -
      0.4*tracegDAdjgD*Sqr(g1)*TSigmaL + 1.2*tracehEAdjhE*Sqr(g1)*TSigmaL + 0.9
      *tracegDAdjgD*Sqr(g1p)*TSigmaL + 0.3*tracehEAdjhE*Sqr(g1p)*TSigmaL + 0.72
      *Sqr(g1)*Sqr(g1p)*TSigmaL + 1.8*Sqr(g1)*Sqr(g2)*TSigmaL + 1.2*Sqr(g1p)*
      Sqr(g2)*TSigmaL + 16*tracegDAdjgD*Sqr(g3)*TSigmaL + 0.02*Power(g1p,4)*Sqr
      (QS)*TSigmaL - 8*Sqr(Conj(KappaPr))*Sqr(KappaPr)*TSigmaL - 50*Sqr(Conj(
      SigmaL))*Sqr(SigmaL)*TSigmaL - 2*Sigmax*Sqr(Conj(Sigmax))*(4*SigmaL*
      TSigmax + Sigmax*TSigmaL) - 0.2*AbsSqr(SigmaL)*(2*SigmaL*(45*
      traceAdjgDTgD + 15*traceAdjhEThE + 6*MassB*Sqr(g1) + 4*MassBp*Sqr(g1p) +
      30*MassWB*Sqr(g2)) - 3*(-45*tracegDAdjgD - 15*tracehEAdjhE + 6*Sqr(g1) +
      4*Sqr(g1p) + 30*Sqr(g2))*TSigmaL + 60*Conj(KappaPr)*(2*SigmaL*TKappaPr +
      3*KappaPr*TSigmaL)) - 0.1*Conj(Sigmax)*(60*traceAdjKappaTKappa*Sigmax*
      SigmaL + 40*traceAdjLambda12TLambda12*Sigmax*SigmaL + 2*MassBp*Sigmax*
      SigmaL*Sqr(g1p)*Sqr(QS) + 60*traceKappaAdjKappa*SigmaL*TSigmax + 40*
      traceLambda12AdjLambda12*SigmaL*TSigmax - 2*SigmaL*Sqr(g1p)*Sqr(QS)*
      TSigmax + 40*Conj(SigmaL)*Sqr(SigmaL)*TSigmax + 30*traceKappaAdjKappa*
      Sigmax*TSigmaL + 20*traceLambda12AdjLambda12*Sigmax*TSigmaL + 60*AbsSqr(
      SigmaL)*Sigmax*TSigmaL - Sigmax*Sqr(g1p)*Sqr(QS)*TSigmaL + 40*Conj(
      KappaPr)*(2*Sigmax*SigmaL*TKappaPr + 2*KappaPr*SigmaL*TSigmax + KappaPr*
      Sigmax*TSigmaL) + 20*Conj(Lambdax)*(2*Sigmax*SigmaL*TLambdax + 2*Lambdax*
      SigmaL*TSigmax + Lambdax*Sigmax*TSigmaL))));


   return beta_TSigmaL;
}

} // namespace flexiblesusy
