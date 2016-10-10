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

// File generated at Wed 3 Jun 2015 23:43:03

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TSigmax.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TSigmax_one_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_TSigmax;

   beta_TSigmax = Re(oneOver16PiSqr*(6*traceAdjKappaTKappa*Sigmax + 4*
      traceAdjLambda12TLambda12*Sigmax + 0.2*MassBp*Sigmax*Sqr(g1p)*Sqr(QS) + 3
      *traceKappaAdjKappa*TSigmax + 2*traceLambda12AdjLambda12*TSigmax + 9*
      AbsSqr(Sigmax)*TSigmax + 2*AbsSqr(SigmaL)*TSigmax - 0.1*Sqr(g1p)*Sqr(QS)*
      TSigmax + 2*Conj(KappaPr)*(2*Sigmax*TKappaPr + KappaPr*TSigmax) + 2*Conj(
      Lambdax)*(2*Sigmax*TLambdax + Lambdax*TSigmax) + 4*Conj(SigmaL)*Sigmax*
      TSigmaL));


   return beta_TSigmax;
}

/**
 * Calculates the two-loop beta function of TSigmax.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TSigmax_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracefuAdjLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12Adjfu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE;
   const double traceKappaAdjgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjfuTfuAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjfuTfuAdjLambda12;
   const double traceLambda12AdjhEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjLambda12;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;
   const double traceAdjfdTfdconjLambda12TpLambda12 =
      TRACE_STRUCT.traceAdjfdTfdconjLambda12TpLambda12;
   const double traceAdjLambda12TpfdconjfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjfdTLambda12;


   double beta_TSigmax;

   beta_TSigmax = Re(twoLoop*(-0.04*Power(g1p,4)*MassBp*Power(QS,4)*
      Sigmax - 4*traceAdjfdTfdconjLambda12TpLambda12*Sigmax - 4*
      traceAdjLambda12TpfdconjfdTLambda12*Sigmax - 4*
      tracefuAdjLambda12TLambda12Adjfu*Sigmax - 12*tracegDAdjKappaTKappaAdjgD*
      Sigmax - 4*tracehEAdjLambda12TLambda12AdjhE*Sigmax - 12*
      traceKappaAdjgDTgDAdjKappa*Sigmax - 24*traceKappaAdjKappaTKappaAdjKappa*
      Sigmax - 4*traceLambda12AdjfuTfuAdjLambda12*Sigmax - 4*
      traceLambda12AdjhEThEAdjLambda12*Sigmax - 16*
      traceLambda12AdjLambda12TLambda12AdjLambda12*Sigmax - 12*traceAdjgDTgD*
      AbsSqr(SigmaL)*Sigmax - 4*traceAdjhEThE*AbsSqr(SigmaL)*Sigmax + 1.6*
      traceAdjKappaTKappa*Sigmax*Sqr(g1) + 2.4*traceAdjLambda12TLambda12*Sigmax
      *Sqr(g1) - 1.6*MassB*traceKappaAdjKappa*Sigmax*Sqr(g1) - 2.4*MassB*
      traceLambda12AdjLambda12*Sigmax*Sqr(g1) - 2.4*MassB*AbsSqr(SigmaL)*Sigmax
      *Sqr(g1) + 3.9*traceAdjKappaTKappa*Sigmax*Sqr(g1p) + 2.6*
      traceAdjLambda12TLambda12*Sigmax*Sqr(g1p) - 3.9*MassBp*traceKappaAdjKappa
      *Sigmax*Sqr(g1p) - 2.6*MassBp*traceLambda12AdjLambda12*Sigmax*Sqr(g1p) -
      1.6*MassBp*AbsSqr(SigmaL)*Sigmax*Sqr(g1p) + 12*traceAdjLambda12TLambda12*
      Sigmax*Sqr(g2) - 12*MassWB*traceLambda12AdjLambda12*Sigmax*Sqr(g2) - 12*
      MassWB*AbsSqr(SigmaL)*Sigmax*Sqr(g2) + 32*traceAdjKappaTKappa*Sigmax*Sqr(
      g3) - 32*MassG*traceKappaAdjKappa*Sigmax*Sqr(g3) - 3.76*Power(g1p,4)*
      MassBp*Sigmax*Sqr(QS) - 0.3*traceAdjKappaTKappa*Sigmax*Sqr(g1p)*Sqr(QS) -
      0.2*traceAdjLambda12TLambda12*Sigmax*Sqr(g1p)*Sqr(QS) + 0.3*MassBp*
      traceKappaAdjKappa*Sigmax*Sqr(g1p)*Sqr(QS) + 0.2*MassBp*
      traceLambda12AdjLambda12*Sigmax*Sqr(g1p)*Sqr(QS) - 16*AbsSqr(SigmaL)*Conj
      (KappaPr)*Sigmax*TKappaPr - 32*KappaPr*Sigmax*Sqr(Conj(KappaPr))*TKappaPr
      + 0.01*Power(g1p,4)*Power(QS,4)*TSigmax - 2*
      tracefuAdjLambda12Lambda12Adjfu*TSigmax - 6*tracegDAdjKappaKappaAdjgD*
      TSigmax - 2*tracehEAdjLambda12Lambda12AdjhE*TSigmax - 6*
      traceKappaAdjKappaKappaAdjKappa*TSigmax - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TSigmax - 2*
      traceLambda12AdjLambda12Tpfdconjfd*TSigmax - 6*tracegDAdjgD*AbsSqr(SigmaL
      )*TSigmax - 2*tracehEAdjhE*AbsSqr(SigmaL)*TSigmax - 8*AbsSqr(KappaPr)*
      AbsSqr(SigmaL)*TSigmax + 0.8*traceKappaAdjKappa*Sqr(g1)*TSigmax + 1.2*
      traceLambda12AdjLambda12*Sqr(g1)*TSigmax + 1.2*AbsSqr(SigmaL)*Sqr(g1)*
      TSigmax + 1.95*traceKappaAdjKappa*Sqr(g1p)*TSigmax + 1.3*
      traceLambda12AdjLambda12*Sqr(g1p)*TSigmax + 0.8*AbsSqr(SigmaL)*Sqr(g1p)*
      TSigmax + 6*traceLambda12AdjLambda12*Sqr(g2)*TSigmax + 6*AbsSqr(SigmaL)*
      Sqr(g2)*TSigmax + 16*traceKappaAdjKappa*Sqr(g3)*TSigmax + 0.94*Power(g1p,
      4)*Sqr(QS)*TSigmax - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*TSigmax -
      0.1*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS)*TSigmax - 8*Sqr(Conj(
      KappaPr))*Sqr(KappaPr)*TSigmax - 30*Sqr(Conj(Sigmax))*Sqr(Sigmax)*TSigmax
      - 4*Sqr(Conj(SigmaL))*Sqr(SigmaL)*TSigmax - 4*Lambdax*Sqr(Conj(Lambdax))
      *(4*Sigmax*TLambdax + Lambdax*TSigmax) - 0.1*Conj(Lambdax)*(2*Sigmax*(20*
      tracefdAdjfd + 20*tracefuAdjfu + 60*traceYdAdjYd + 20*traceYeAdjYe + 60*
      traceYuAdjYu + 40*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2)
      + Sqr(g1p)*Sqr(QS))*TLambdax + Lambdax*(2*Sigmax*(20*traceAdjfdTfd + 20*
      traceAdjfuTfu + 60*traceAdjYdTYd + 20*traceAdjYeTYe + 60*traceAdjYuTYu +
      12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr(g2) - MassBp*Sqr(
      g1p)*Sqr(QS)) + (20*tracefdAdjfd + 20*tracefuAdjfu + 60*traceYdAdjYd + 20
      *traceYeAdjYe + 60*traceYuAdjYu + 120*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*
      Sqr(g1p) - 60*Sqr(g2) + Sqr(g1p)*Sqr(QS))*TSigmax)) - 12*tracegDAdjgD*
      Conj(SigmaL)*Sigmax*TSigmaL - 4*tracehEAdjhE*Conj(SigmaL)*Sigmax*TSigmaL
      - 16*AbsSqr(KappaPr)*Conj(SigmaL)*Sigmax*TSigmaL + 2.4*Conj(SigmaL)*
      Sigmax*Sqr(g1)*TSigmaL + 1.6*Conj(SigmaL)*Sigmax*Sqr(g1p)*TSigmaL + 12*
      Conj(SigmaL)*Sigmax*Sqr(g2)*TSigmaL - 16*Sigmax*SigmaL*Sqr(Conj(SigmaL))*
      TSigmaL - 0.1*AbsSqr(Sigmax)*(-3*(-60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 40*AbsSqr(SigmaL) + Sqr(g1p)*Sqr(QS))*TSigmax
      + 80*Conj(KappaPr)*(2*Sigmax*TKappaPr + 3*KappaPr*TSigmax) + 2*Sigmax*(60
      *traceAdjKappaTKappa + 40*traceAdjLambda12TLambda12 + MassBp*Sqr(g1p)*Sqr
      (QS) + 40*Conj(SigmaL)*TSigmaL))));


   return beta_TSigmax;
}

} // namespace flexiblesusy
