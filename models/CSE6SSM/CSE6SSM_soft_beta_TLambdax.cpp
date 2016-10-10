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

// File generated at Wed 3 Jun 2015 23:43:10

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*((tracefdAdjfd + tracefuAdjfu + 3*
      traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) + AbsSqr(Sigmax) - 0.6
      *Sqr(g1) - 0.65*Sqr(g1p) - 3*Sqr(g2) - 0.05*Sqr(g1p)*Sqr(QS))*TLambdax +
      0.1*Lambdax*(20*traceAdjfdTfd + 20*traceAdjfuTfu + 60*traceAdjKappaTKappa
      + 40*traceAdjLambda12TLambda12 + 60*traceAdjYdTYd + 20*traceAdjYeTYe +
      60*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr(
      g2) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax)));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double tracefdAdjfdfdAdjfd = TRACE_STRUCT.tracefdAdjfdfdAdjfd;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefdAdjfdTfdAdjfd = TRACE_STRUCT.tracefdAdjfdTfdAdjfd;
   const double tracefdAdjfdTfuAdjfu = TRACE_STRUCT.tracefdAdjfdTfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjfuTfdAdjfd = TRACE_STRUCT.tracefuAdjfuTfdAdjfd;
   const double tracefuAdjfuTfuAdjfu = TRACE_STRUCT.tracefuAdjfuTfuAdjfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjhEThEAdjfu = TRACE_STRUCT.tracefuAdjhEThEAdjfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracefuAdjLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12Adjfu;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjfuTfuAdjhE = TRACE_STRUCT.tracehEAdjfuTfuAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjhETYeAdjYe = TRACE_STRUCT.tracehEAdjhETYeAdjYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeThEAdjhE = TRACE_STRUCT.traceYeAdjYeThEAdjhE;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
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
   const double traceAdjgDTpYdconjYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjYdTgD;
   const double traceAdjgDTpYuconjYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjYuTgD;
   const double traceAdjYdTYdconjgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjgDTpgD;
   const double traceAdjYuTYuconjgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjgDTpgD;
   const double traceAdjLambda12TpfdconjfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjfdTLambda12;


   double beta_TLambdax;

   beta_TLambdax = Re(twoLoop*(-23.76*Power(g1,4)*MassB*Lambdax - 25.41*
      Power(g1p,4)*MassBp*Lambdax - 66*Power(g2,4)*MassWB*Lambdax - 0.02*Power(
      g1p,4)*MassBp*Power(QS,4)*Lambdax - 6*traceAdjfdTfdconjLambda12TpLambda12
      *Lambdax - 6*traceAdjgDTpYdconjYdTgD*Lambdax - 6*traceAdjgDTpYuconjYuTgD*
      Lambdax - 6*traceAdjLambda12TpfdconjfdTLambda12*Lambdax - 6*
      traceAdjYdTYdconjgDTpgD*Lambdax - 6*traceAdjYuTYuconjgDTpgD*Lambdax - 12*
      tracefdAdjfdTfdAdjfd*Lambdax - 8*tracefdAdjfdTfuAdjfu*Lambdax - 8*
      tracefuAdjfuTfdAdjfd*Lambdax - 12*tracefuAdjfuTfuAdjfu*Lambdax - 2*
      tracefuAdjhEThEAdjfu*Lambdax - 6*tracefuAdjLambda12TLambda12Adjfu*Lambdax
      - 12*tracegDAdjKappaTKappaAdjgD*Lambdax - 2*tracehEAdjfuTfuAdjhE*Lambdax
      - 4*tracehEAdjhETYeAdjYe*Lambdax - 4*tracehEAdjLambda12TLambda12AdjhE*
      Lambdax - 12*traceKappaAdjgDTgDAdjKappa*Lambdax - 24*
      traceKappaAdjKappaTKappaAdjKappa*Lambdax - 6*
      traceLambda12AdjfuTfuAdjLambda12*Lambdax - 4*
      traceLambda12AdjhEThEAdjLambda12*Lambdax - 16*
      traceLambda12AdjLambda12TLambda12AdjLambda12*Lambdax - 36*
      traceYdAdjYdTYdAdjYd*Lambdax - 12*traceYdAdjYuTYuAdjYd*Lambdax - 4*
      traceYeAdjYeThEAdjhE*Lambdax - 12*traceYeAdjYeTYeAdjYe*Lambdax - 12*
      traceYuAdjYdTYdAdjYu*Lambdax - 36*traceYuAdjYuTYuAdjYu*Lambdax + 1.6*
      traceAdjKappaTKappa*Lambdax*Sqr(g1) + 2.4*traceAdjLambda12TLambda12*
      Lambdax*Sqr(g1) - 0.8*traceAdjYdTYd*Lambdax*Sqr(g1) + 2.4*traceAdjYeTYe*
      Lambdax*Sqr(g1) + 1.6*traceAdjYuTYu*Lambdax*Sqr(g1) - 1.6*MassB*
      traceKappaAdjKappa*Lambdax*Sqr(g1) - 2.4*MassB*traceLambda12AdjLambda12*
      Lambdax*Sqr(g1) + 0.8*MassB*traceYdAdjYd*Lambdax*Sqr(g1) - 2.4*MassB*
      traceYeAdjYe*Lambdax*Sqr(g1) - 1.6*MassB*traceYuAdjYu*Lambdax*Sqr(g1) + 2
      *traceAdjfdTfd*Lambdax*Sqr(g1p) + 3*traceAdjfuTfu*Lambdax*Sqr(g1p) + 3.9*
      traceAdjKappaTKappa*Lambdax*Sqr(g1p) + 2.6*traceAdjLambda12TLambda12*
      Lambdax*Sqr(g1p) - 1.2*traceAdjYdTYd*Lambdax*Sqr(g1p) - 0.4*traceAdjYeTYe
      *Lambdax*Sqr(g1p) - 0.6*traceAdjYuTYu*Lambdax*Sqr(g1p) - 2*MassBp*
      tracefdAdjfd*Lambdax*Sqr(g1p) - 3*MassBp*tracefuAdjfu*Lambdax*Sqr(g1p) -
      3.9*MassBp*traceKappaAdjKappa*Lambdax*Sqr(g1p) - 2.6*MassBp*
      traceLambda12AdjLambda12*Lambdax*Sqr(g1p) + 1.2*MassBp*traceYdAdjYd*
      Lambdax*Sqr(g1p) + 0.4*MassBp*traceYeAdjYe*Lambdax*Sqr(g1p) + 0.6*MassBp*
      traceYuAdjYu*Lambdax*Sqr(g1p) - 0.54*MassB*Lambdax*Sqr(g1)*Sqr(g1p) -
      0.54*MassBp*Lambdax*Sqr(g1)*Sqr(g1p) + 12*traceAdjLambda12TLambda12*
      Lambdax*Sqr(g2) - 12*MassWB*traceLambda12AdjLambda12*Lambdax*Sqr(g2) -
      3.6*MassB*Lambdax*Sqr(g1)*Sqr(g2) - 3.6*MassWB*Lambdax*Sqr(g1)*Sqr(g2) -
      3.9*MassBp*Lambdax*Sqr(g1p)*Sqr(g2) - 3.9*MassWB*Lambdax*Sqr(g1p)*Sqr(g2)
      + 32*traceAdjKappaTKappa*Lambdax*Sqr(g3) + 32*traceAdjYdTYd*Lambdax*Sqr(
      g3) + 32*traceAdjYuTYu*Lambdax*Sqr(g3) - 32*MassG*traceKappaAdjKappa*
      Lambdax*Sqr(g3) - 32*MassG*traceYdAdjYd*Lambdax*Sqr(g3) - 32*MassG*
      traceYuAdjYu*Lambdax*Sqr(g3) - 2.01*Power(g1p,4)*MassBp*Lambdax*Sqr(QS) -
      0.3*traceAdjKappaTKappa*Lambdax*Sqr(g1p)*Sqr(QS) - 0.2*
      traceAdjLambda12TLambda12*Lambdax*Sqr(g1p)*Sqr(QS) + 0.3*MassBp*
      traceKappaAdjKappa*Lambdax*Sqr(g1p)*Sqr(QS) + 0.2*MassBp*
      traceLambda12AdjLambda12*Lambdax*Sqr(g1p)*Sqr(QS) + 5.94*Power(g1,4)*
      TLambdax + 6.3525*Power(g1p,4)*TLambdax + 16.5*Power(g2,4)*TLambdax +
      0.005*Power(g1p,4)*Power(QS,4)*TLambdax - 3*tracefdAdjfdfdAdjfd*TLambdax
      - 4*tracefdAdjfdfuAdjfu*TLambdax - 3*tracefuAdjfufuAdjfu*TLambdax -
      tracefuAdjhEhEAdjfu*TLambdax - 3*tracefuAdjLambda12Lambda12Adjfu*TLambdax
      - 3*tracegDAdjgDTpYdconjYd*TLambdax - 3*tracegDAdjgDTpYuconjYu*TLambdax
      - 6*tracegDAdjKappaKappaAdjgD*TLambdax - 2*tracehEAdjhEYeAdjYe*TLambdax -
      2*tracehEAdjLambda12Lambda12AdjhE*TLambdax - 6*
      traceKappaAdjKappaKappaAdjKappa*TLambdax - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TLambdax - 3*
      traceLambda12AdjLambda12Tpfdconjfd*TLambdax - 9*traceYdAdjYdYdAdjYd*
      TLambdax - 6*traceYdAdjYuYuAdjYd*TLambdax - 3*traceYeAdjYeYeAdjYe*
      TLambdax - 9*traceYuAdjYuYuAdjYu*TLambdax - 2*AbsSqr(Sigmax)*AbsSqr(
      SigmaL)*TLambdax + 0.8*traceKappaAdjKappa*Sqr(g1)*TLambdax + 1.2*
      traceLambda12AdjLambda12*Sqr(g1)*TLambdax - 0.4*traceYdAdjYd*Sqr(g1)*
      TLambdax + 1.2*traceYeAdjYe*Sqr(g1)*TLambdax + 0.8*traceYuAdjYu*Sqr(g1)*
      TLambdax + tracefdAdjfd*Sqr(g1p)*TLambdax + 1.5*tracefuAdjfu*Sqr(g1p)*
      TLambdax + 1.95*traceKappaAdjKappa*Sqr(g1p)*TLambdax + 1.3*
      traceLambda12AdjLambda12*Sqr(g1p)*TLambdax - 0.6*traceYdAdjYd*Sqr(g1p)*
      TLambdax - 0.2*traceYeAdjYe*Sqr(g1p)*TLambdax - 0.3*traceYuAdjYu*Sqr(g1p)
      *TLambdax + 0.27*Sqr(g1)*Sqr(g1p)*TLambdax + 6*traceLambda12AdjLambda12*
      Sqr(g2)*TLambdax + 1.8*Sqr(g1)*Sqr(g2)*TLambdax + 1.95*Sqr(g1p)*Sqr(g2)*
      TLambdax + 16*traceKappaAdjKappa*Sqr(g3)*TLambdax + 16*traceYdAdjYd*Sqr(
      g3)*TLambdax + 16*traceYuAdjYu*Sqr(g3)*TLambdax + 0.5025*Power(g1p,4)*Sqr
      (QS)*TLambdax - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*TLambdax - 0.1*
      traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS)*TLambdax - 50*Sqr(Conj(Lambdax)
      )*Sqr(Lambdax)*TLambdax - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)*TLambdax - 4*
      AbsSqr(SigmaL)*Conj(Sigmax)*Lambdax*TSigmax - 8*Lambdax*Sigmax*Sqr(Conj(
      Sigmax))*TSigmax - 2*Conj(KappaPr)*Conj(Sigmax)*(2*Lambdax*Sigmax*
      TKappaPr + KappaPr*Sigmax*TLambdax + 2*KappaPr*Lambdax*TSigmax) - 0.1*
      AbsSqr(Lambdax)*(3*(30*tracefdAdjfd + 30*tracefuAdjfu + 60*
      traceKappaAdjKappa + 40*traceLambda12AdjLambda12 + 90*traceYdAdjYd + 30*
      traceYeAdjYe + 90*traceYuAdjYu + 20*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(
      g1p) - 60*Sqr(g2))*TLambdax + 2*Lambdax*(30*traceAdjfdTfd + 30*
      traceAdjfuTfu + 60*traceAdjKappaTKappa + 40*traceAdjLambda12TLambda12 +
      90*traceAdjYdTYd + 30*traceAdjYeTYe + 90*traceAdjYuTYu + 12*MassB*Sqr(g1)
      + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr(g2) + 20*Conj(Sigmax)*TSigmax)) - 4
      *AbsSqr(Sigmax)*Conj(SigmaL)*Lambdax*TSigmaL));


   return beta_TLambdax;
}

} // namespace flexiblesusy
