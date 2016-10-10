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

// File generated at Wed 3 Jun 2015 23:43:25

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTfuTpTfu = TRACE_STRUCT.traceconjTfuTpTfu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracefumH1I2Adjfu = TRACE_STRUCT.tracefumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2 = TRACE_STRUCT.tracefuAdjfuconjmSI2;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 -
      0.6324555320336759*g1p*Tr14 + 2*traceconjTfuTpTfu + 6*traceconjTYuTpTYu +
      2*mHu2*tracefuAdjfu + 2*tracefuAdjfuconjmSI2 + 2*tracefumH1I2Adjfu + 6*
      tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*traceYuAdjYu + 2*mHd2*AbsSqr
      (Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(
      TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 0.8*AbsSqr(MassBp)*Sqr(g1p) - 6*
      AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTfdTpfd = TRACE_STRUCT.traceconjTfdTpfd;
   const double traceconjTfdTpTfd = TRACE_STRUCT.traceconjTfdTpTfd;
   const double traceconjTfuTpfu = TRACE_STRUCT.traceconjTfuTpfu;
   const double traceconjTfuTpTfu = TRACE_STRUCT.traceconjTfuTpTfu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracefdmH2I2Adjfd = TRACE_STRUCT.tracefdmH2I2Adjfd;
   const double tracefdAdjfdconjmSI2 = TRACE_STRUCT.tracefdAdjfdconjmSI2;
   const double tracefumH1I2Adjfu = TRACE_STRUCT.tracefumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2 = TRACE_STRUCT.tracefuAdjfuconjmSI2;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefdAdjfdTfuAdjTfu =
      TRACE_STRUCT.tracefdAdjfdTfuAdjTfu;
   const double tracefdAdjTfdTfuAdjfu =
      TRACE_STRUCT.tracefdAdjTfdTfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjfuTfdAdjTfd =
      TRACE_STRUCT.tracefuAdjfuTfdAdjTfd;
   const double tracefuAdjfuTfuAdjTfu =
      TRACE_STRUCT.tracefuAdjfuTfuAdjTfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjhEThEAdjTfu =
      TRACE_STRUCT.tracefuAdjhEThEAdjTfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracefuAdjLambda12TLambda12AdjTfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12AdjTfu;
   const double tracefuAdjTfuTfdAdjfd =
      TRACE_STRUCT.tracefuAdjTfuTfdAdjfd;
   const double tracefuAdjTfuTfuAdjfu =
      TRACE_STRUCT.tracefuAdjTfuTfuAdjfu;
   const double tracefuAdjThEThEAdjfu =
      TRACE_STRUCT.tracefuAdjThEThEAdjfu;
   const double tracefuAdjTLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjTLambda12TLambda12Adjfu;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjgDTpTYuconjTYu =
      TRACE_STRUCT.tracegDAdjgDTpTYuconjTYu;
   const double tracehEAdjfuTfuAdjThE =
      TRACE_STRUCT.tracehEAdjfuTfuAdjThE;
   const double tracehEAdjTfuTfuAdjhE =
      TRACE_STRUCT.tracehEAdjTfuTfuAdjhE;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double traceYuconjTgDTpTgDAdjYu =
      TRACE_STRUCT.traceYuconjTgDTpTgDAdjYu;
   const double traceLambda12AdjfuTfuAdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjfuTfuAdjTLambda12;
   const double traceLambda12AdjTfuTfuAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTfuTfuAdjLambda12;
   const double traceAdjgDTpYuconjTYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjTYuTgD;
   const double traceAdjYuTYuconjTgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjTgDTpgD;
   const double tracefdmH2I2AdjfdfuAdjfu =
      TRACE_STRUCT.tracefdmH2I2AdjfdfuAdjfu;
   const double tracefdAdjfdfumH1I2Adjfu =
      TRACE_STRUCT.tracefdAdjfdfumH1I2Adjfu;
   const double tracefdAdjfdfuAdjfuconjmSI2 =
      TRACE_STRUCT.tracefdAdjfdfuAdjfuconjmSI2;
   const double tracefdAdjfdconjmSI2fuAdjfu =
      TRACE_STRUCT.tracefdAdjfdconjmSI2fuAdjfu;
   const double tracefumH1I2AdjfufuAdjfu =
      TRACE_STRUCT.tracefumH1I2AdjfufuAdjfu;
   const double tracefumH1I2AdjhEhEAdjfu =
      TRACE_STRUCT.tracefumH1I2AdjhEhEAdjfu;
   const double tracefumH1I2AdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefumH1I2AdjLambda12Lambda12Adjfu;
   const double tracefuAdjfufumH1I2Adjfu =
      TRACE_STRUCT.tracefuAdjfufumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2fuAdjfu =
      TRACE_STRUCT.tracefuAdjfuconjmSI2fuAdjfu;
   const double tracefuAdjhEhEmH1I2Adjfu =
      TRACE_STRUCT.tracefuAdjhEhEmH1I2Adjfu;
   const double tracefuAdjhEhEAdjfuconjmSI2 =
      TRACE_STRUCT.tracefuAdjhEhEAdjfuconjmSI2;
   const double tracefuAdjhEme2hEAdjfu =
      TRACE_STRUCT.tracefuAdjhEme2hEAdjfu;
   const double tracefuAdjLambda12Lambda12mH1I2Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12mH1I2Adjfu;
   const double tracefuAdjLambda12Lambda12AdjfuconjmSI2 =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12AdjfuconjmSI2;
   const double tracefuAdjLambda12conjmH2I2Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12conjmH2I2Lambda12Adjfu;
   const double tracegDAdjgDconjmq2TpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDconjmq2TpYuconjYu;
   const double tracegDAdjgDTpYuconjmu2conjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjmu2conjYu;
   const double tracegDAdjgDTpYuconjYuconjmq2 =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYuconjmq2;
   const double tracegDconjmDxbar2AdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYuconjYu;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHu2;

   beta_mHu2 = Re(twoLoop*(6*Power(g2,4)*Tr22 - 0.9797958971132712*g1*g1p
      *Tr2U114 - 0.9797958971132712*g1*g1p*Tr2U141 + 3.0983866769659336*g1*Tr31
      - 2.5298221281347035*g1p*Tr34 - 6*traceAdjgDTpYuconjTYuTgD - 6*
      traceAdjYuTYuconjTgDTpgD - 4*tracefdAdjfdconjmSI2fuAdjfu - 4*mHd2*
      tracefdAdjfdfuAdjfu - 4*mHu2*tracefdAdjfdfuAdjfu - 4*
      tracefdAdjfdfuAdjfuconjmSI2 - 4*tracefdAdjfdfumH1I2Adjfu - 4*
      tracefdAdjfdTfuAdjTfu - 4*tracefdAdjTfdTfuAdjfu - 4*
      tracefdmH2I2AdjfdfuAdjfu - 12*tracefuAdjfuconjmSI2fuAdjfu - 12*mHu2*
      tracefuAdjfufuAdjfu - 6*tracefuAdjfufumH1I2Adjfu - 4*
      tracefuAdjfuTfdAdjTfd - 12*tracefuAdjfuTfuAdjTfu - 2*mHp2*
      tracefuAdjhEhEAdjfu - 2*mHu2*tracefuAdjhEhEAdjfu - 2*
      tracefuAdjhEhEAdjfuconjmSI2 - 2*tracefuAdjhEhEmH1I2Adjfu - 2*
      tracefuAdjhEme2hEAdjfu - 2*tracefuAdjhEThEAdjTfu - 2*
      tracefuAdjLambda12conjmH2I2Lambda12Adjfu - 2*mHu2*
      tracefuAdjLambda12Lambda12Adjfu - 2*ms2*tracefuAdjLambda12Lambda12Adjfu -
      2*tracefuAdjLambda12Lambda12AdjfuconjmSI2 - 2*
      tracefuAdjLambda12Lambda12mH1I2Adjfu - 2*
      tracefuAdjLambda12TLambda12AdjTfu - 4*tracefuAdjTfuTfdAdjfd - 12*
      tracefuAdjTfuTfuAdjfu - 2*tracefuAdjThEThEAdjfu - 2*
      tracefuAdjTLambda12TLambda12Adjfu - 6*tracefumH1I2AdjfufuAdjfu - 2*
      tracefumH1I2AdjhEhEAdjfu - 2*tracefumH1I2AdjLambda12Lambda12Adjfu - 6*
      tracegDAdjgDconjmq2TpYuconjYu - 6*tracegDAdjgDTpTYuconjTYu - 6*
      tracegDAdjgDTpYuconjmu2conjYu - 6*mHp2*tracegDAdjgDTpYuconjYu - 6*mHu2*
      tracegDAdjgDTpYuconjYu - 6*tracegDAdjgDTpYuconjYuconjmq2 - 6*
      tracegDconjmDxbar2AdjgDTpYuconjYu - 2*tracehEAdjfuTfuAdjThE - 2*
      tracehEAdjTfuTfuAdjhE - 2*traceLambda12AdjfuTfuAdjTLambda12 - 2*
      traceLambda12AdjTfuTfuAdjLambda12 - 6*tracemd2YdAdjYuYuAdjYd - 6*
      tracemq2AdjYdYdAdjYuYu - 6*tracemq2AdjYuYuAdjYdYd - 36*
      tracemq2AdjYuYuAdjYuYu - 6*tracemu2YuAdjYdYdAdjYu - 36*
      tracemu2YuAdjYuYuAdjYu - 6*traceYdAdjTYuTYuAdjYd - 6*
      traceYdAdjYuTYuAdjTYd - 6*mHd2*traceYdAdjYuYuAdjYd - 6*mHu2*
      traceYdAdjYuYuAdjYd - 6*traceYuAdjTYdTYdAdjYu - 36*traceYuAdjTYuTYuAdjYu
      - 6*traceYuAdjYdTYdAdjTYu - 36*traceYuAdjYuTYuAdjTYu - 36*mHu2*
      traceYuAdjYuYuAdjYu - 6*traceYuconjTgDTpTgDAdjYu + 87*Power(g2,4)*AbsSqr(
      MassWB) - 2*traceconjTfdTpTfd*AbsSqr(Lambdax) - 6*traceconjTKappaTpTKappa
      *AbsSqr(Lambdax) - 4*traceconjTLambda12TpTLambda12*AbsSqr(Lambdax) - 6*
      traceconjTYdTpTYd*AbsSqr(Lambdax) - 2*traceconjTYeTpTYe*AbsSqr(Lambdax) -
      4*mHd2*tracefdAdjfd*AbsSqr(Lambdax) - 2*mHu2*tracefdAdjfd*AbsSqr(Lambdax
      ) - 2*ms2*tracefdAdjfd*AbsSqr(Lambdax) - 2*tracefdAdjfdconjmSI2*AbsSqr(
      Lambdax) - 2*tracefdmH2I2Adjfd*AbsSqr(Lambdax) - 6*mHd2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 6*mHu2*traceKappaAdjKappa*AbsSqr(
      Lambdax) - 12*ms2*traceKappaAdjKappa*AbsSqr(Lambdax) - 6*
      traceKappaAdjKappaconjmDx2*AbsSqr(Lambdax) - 6*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 4*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 4*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 8*ms2*traceLambda12AdjLambda12
      *AbsSqr(Lambdax) - 4*traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) -
      6*tracemd2YdAdjYd*AbsSqr(Lambdax) - 2*traceme2YeAdjYe*AbsSqr(Lambdax) - 4
      *tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 2*traceml2AdjYeYe*AbsSqr
      (Lambdax) - 6*tracemq2AdjYdYd*AbsSqr(Lambdax) - 12*mHd2*traceYdAdjYd*
      AbsSqr(Lambdax) - 6*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 6*ms2*
      traceYdAdjYd*AbsSqr(Lambdax) - 4*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 2*
      mHu2*traceYeAdjYe*AbsSqr(Lambdax) - 2*ms2*traceYeAdjYe*AbsSqr(Lambdax) -
      2*mHd2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*mHu2*AbsSqr(Lambdax)*AbsSqr(
      Sigmax) - 2*mphi2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 4*ms2*AbsSqr(Lambdax)*
      AbsSqr(Sigmax) - 2*msbar2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*tracefdAdjfd
      *AbsSqr(TLambdax) - 6*traceKappaAdjKappa*AbsSqr(TLambdax) - 4*
      traceLambda12AdjLambda12*AbsSqr(TLambdax) - 6*traceYdAdjYd*AbsSqr(
      TLambdax) - 2*traceYeAdjYe*AbsSqr(TLambdax) - 24*AbsSqr(Lambdax)*AbsSqr(
      TLambdax) - 2*AbsSqr(Sigmax)*AbsSqr(TLambdax) - 2*AbsSqr(Lambdax)*AbsSqr(
      TSigmax) - 2*traceAdjfdTfd*Conj(TLambdax)*Lambdax - 6*traceAdjKappaTKappa
      *Conj(TLambdax)*Lambdax - 4*traceAdjLambda12TLambda12*Conj(TLambdax)*
      Lambdax - 6*traceAdjYdTYd*Conj(TLambdax)*Lambdax - 2*traceAdjYeTYe*Conj(
      TLambdax)*Lambdax + 1.2*Tr2U111*Sqr(g1) + 1.6*traceconjTYuTpTYu*Sqr(g1) -
      1.6*MassB*traceconjTYuTpYu*Sqr(g1) + 1.6*tracemq2AdjYuYu*Sqr(g1) + 1.6*
      tracemu2YuAdjYu*Sqr(g1) + 1.6*mHu2*traceYuAdjYu*Sqr(g1) + 0.8*Tr2U144*Sqr
      (g1p) - 3*MassBp*traceconjTfuTpfu*Sqr(g1p) + 3*traceconjTfuTpTfu*Sqr(g1p)
      - 0.6*traceconjTYuTpTYu*Sqr(g1p) + 0.6*MassBp*traceconjTYuTpYu*Sqr(g1p)
      + 3*mHu2*tracefuAdjfu*Sqr(g1p) + 3*tracefuAdjfuconjmSI2*Sqr(g1p) + 3*
      tracefumH1I2Adjfu*Sqr(g1p) - 0.6*tracemq2AdjYuYu*Sqr(g1p) - 0.6*
      tracemu2YuAdjYu*Sqr(g1p) - 0.6*mHu2*traceYuAdjYu*Sqr(g1p) + 0.5*mHd2*
      AbsSqr(Lambdax)*Sqr(g1p) + 0.5*mHu2*AbsSqr(Lambdax)*Sqr(g1p) + 0.5*ms2*
      AbsSqr(Lambdax)*Sqr(g1p) + 0.5*AbsSqr(TLambdax)*Sqr(g1p) - 0.5*MassBp*
      Conj(TLambdax)*Lambdax*Sqr(g1p) + 3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) +
      1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 2.4*AbsSqr(MassWB)*Sqr(g1p)*Sqr(
      g2) + 1.2*MassBp*Conj(MassWB)*Sqr(g1p)*Sqr(g2) + 0.04*Conj(MassB)*Sqr(g1)
      *(-40*traceAdjYuTYu + 80*MassB*traceYuAdjYu + 9*(99*MassB*Sqr(g1) + 2*(2*
      MassB + MassBp)*Sqr(g1p) + 5*(2*MassB + MassWB)*Sqr(g2))) + 32*
      traceconjTYuTpTYu*Sqr(g3) - 32*MassG*traceconjTYuTpYu*Sqr(g3) + 32*
      tracemq2AdjYuYu*Sqr(g3) + 32*tracemu2YuAdjYu*Sqr(g3) + 32*mHu2*
      traceYuAdjYu*Sqr(g3) + 64*traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYuTYu*Conj(MassG)*Sqr(g3) + 0.1*mHd2*AbsSqr(Lambdax)*Sqr(g1p)*Sqr
      (QS) + 0.1*mHu2*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS) + 0.1*ms2*AbsSqr(Lambdax
      )*Sqr(g1p)*Sqr(QS) + 0.1*AbsSqr(TLambdax)*Sqr(g1p)*Sqr(QS) - 0.1*MassBp*
      Conj(TLambdax)*Lambdax*Sqr(g1p)*Sqr(QS) - 12*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 12*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*ms2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 0.02*Conj(MassBp)*Sqr(g1p)*(6*(-25*traceAdjfuTfu
      + 5*traceAdjYuTYu + 50*MassBp*tracefuAdjfu - 10*MassBp*traceYuAdjYu + 6*
      MassB*Sqr(g1) + 12*MassBp*Sqr(g1) + 192*MassBp*Sqr(g1p) + 20*MassBp*Sqr(
      g2) + 10*MassWB*Sqr(g2) + MassBp*Sqr(g1p)*Sqr(QS)) + 5*Conj(Lambdax)*(5 +
      Sqr(QS))*(2*MassBp*Lambdax - TLambdax)) - 2*traceconjTfdTpfd*Conj(
      Lambdax)*TLambdax - 6*traceconjTKappaTpKappa*Conj(Lambdax)*TLambdax - 4*
      traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax - 6*traceconjTYdTpYd*
      Conj(Lambdax)*TLambdax - 2*traceconjTYeTpYe*Conj(Lambdax)*TLambdax - 2*
      Conj(Lambdax)*Conj(TSigmax)*Sigmax*TLambdax - 2*Conj(Sigmax)*Conj(
      TLambdax)*Lambdax*TSigmax));


   return beta_mHu2;
}

} // namespace flexiblesusy
