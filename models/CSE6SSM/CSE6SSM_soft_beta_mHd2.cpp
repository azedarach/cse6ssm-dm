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

// File generated at Wed 3 Jun 2015 23:43:22

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTfdTpTfd = TRACE_STRUCT.traceconjTfdTpTfd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracefdmH2I2Adjfd = TRACE_STRUCT.tracefdmH2I2Adjfd;
   const double tracefdAdjfdconjmSI2 = TRACE_STRUCT.tracefdAdjfdconjmSI2;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 -
      0.9486832980505138*g1p*Tr14 + 2*traceconjTfdTpTfd + 6*traceconjTYdTpTYd +
      2*traceconjTYeTpTYe + 2*mHd2*tracefdAdjfd + 2*tracefdAdjfdconjmSI2 + 2*
      tracefdmH2I2Adjfd + 6*tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*
      traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*
      traceYeAdjYe + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*
      AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 1.8*
      AbsSqr(MassBp)*Sqr(g1p) - 6*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
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
   const double tracefdAdjfdfdAdjfd = TRACE_STRUCT.tracefdAdjfdfdAdjfd;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefdAdjfdTfdAdjTfd =
      TRACE_STRUCT.tracefdAdjfdTfdAdjTfd;
   const double tracefdAdjfdTfuAdjTfu =
      TRACE_STRUCT.tracefdAdjfdTfuAdjTfu;
   const double tracefdAdjTfdTfdAdjfd =
      TRACE_STRUCT.tracefdAdjTfdTfdAdjfd;
   const double tracefdAdjTfdTfuAdjfu =
      TRACE_STRUCT.tracefdAdjTfdTfuAdjfu;
   const double tracefdconjTLambda12TpTLambda12Adjfd =
      TRACE_STRUCT.tracefdconjTLambda12TpTLambda12Adjfd;
   const double tracefuAdjfuTfdAdjTfd =
      TRACE_STRUCT.tracefuAdjfuTfdAdjTfd;
   const double tracefuAdjTfuTfdAdjfd =
      TRACE_STRUCT.tracefuAdjTfuTfdAdjfd;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpTYdconjTYd =
      TRACE_STRUCT.tracegDAdjgDTpTYdconjTYd;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjhETYeAdjTYe =
      TRACE_STRUCT.tracehEAdjhETYeAdjTYe;
   const double tracehEAdjThETYeAdjYe =
      TRACE_STRUCT.tracehEAdjThETYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYdconjTgDTpTgDAdjYd =
      TRACE_STRUCT.traceYdconjTgDTpTgDAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeThEAdjThE =
      TRACE_STRUCT.traceYeAdjYeThEAdjThE;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeThEAdjhE =
      TRACE_STRUCT.traceYeAdjTYeThEAdjhE;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;
   const double traceLambda12AdjLambda12TpTfdconjTfd =
      TRACE_STRUCT.traceLambda12AdjLambda12TpTfdconjTfd;
   const double traceAdjfdTfdconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceAdjfdTfdconjTLambda12TpLambda12;
   const double traceAdjgDTpYdconjTYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjTYdTgD;
   const double traceAdjYdTYdconjTgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjTgDTpgD;
   const double traceAdjLambda12TpfdconjTfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjTfdTLambda12;
   const double tracefdmH2I2AdjfdfdAdjfd =
      TRACE_STRUCT.tracefdmH2I2AdjfdfdAdjfd;
   const double tracefdmH2I2AdjfdfuAdjfu =
      TRACE_STRUCT.tracefdmH2I2AdjfdfuAdjfu;
   const double tracefdAdjfdfdmH2I2Adjfd =
      TRACE_STRUCT.tracefdAdjfdfdmH2I2Adjfd;
   const double tracefdAdjfdfumH1I2Adjfu =
      TRACE_STRUCT.tracefdAdjfdfumH1I2Adjfu;
   const double tracefdAdjfdfuAdjfuconjmSI2 =
      TRACE_STRUCT.tracefdAdjfdfuAdjfuconjmSI2;
   const double tracefdAdjfdconjmSI2fdAdjfd =
      TRACE_STRUCT.tracefdAdjfdconjmSI2fdAdjfd;
   const double tracefdAdjfdconjmSI2fuAdjfu =
      TRACE_STRUCT.tracefdAdjfdconjmSI2fuAdjfu;
   const double tracefdconjLambda12TpLambda12AdjfdconjmSI2 =
      TRACE_STRUCT.tracefdconjLambda12TpLambda12AdjfdconjmSI2;
   const double tracegDAdjgDconjmq2TpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDconjmq2TpYdconjYd;
   const double tracegDAdjgDTpYdconjmd2conjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjmd2conjYd;
   const double tracegDAdjgDTpYdconjYdconjmq2 =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYdconjmq2;
   const double tracegDconjmDxbar2AdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYdconjYd;
   const double tracehEmH1I2AdjhEYeAdjYe =
      TRACE_STRUCT.tracehEmH1I2AdjhEYeAdjYe;
   const double tracehEAdjhEme2YeAdjYe =
      TRACE_STRUCT.tracehEAdjhEme2YeAdjYe;
   const double tracehEAdjhEYeml2AdjYe =
      TRACE_STRUCT.tracehEAdjhEYeml2AdjYe;
   const double tracehEAdjhEYeAdjYeme2 =
      TRACE_STRUCT.tracehEAdjhEYeAdjYeme2;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double tracemH1I2AdjLambda12TpfdconjfdLambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12TpfdconjfdLambda12;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double traceLambda12AdjLambda12conjmH2I2Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Tpfdconjfd;
   const double traceLambda12AdjLambda12TpfdconjfdconjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12TpfdconjfdconjmH2I2;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   beta_mHd2 = Re(twoLoop*(6*Power(g2,4)*Tr22 + 1.4696938456699067*g1*g1p
      *Tr2U114 + 1.4696938456699067*g1*g1p*Tr2U141 - 3.0983866769659336*g1*Tr31
      - 3.794733192202055*g1p*Tr34 - 2*traceAdjfdTfdconjTLambda12TpLambda12 -
      6*traceAdjgDTpYdconjTYdTgD - 2*traceAdjLambda12TpfdconjTfdTLambda12 - 6*
      traceAdjYdTYdconjTgDTpgD - 12*tracefdAdjfdconjmSI2fdAdjfd - 4*
      tracefdAdjfdconjmSI2fuAdjfu - 12*mHd2*tracefdAdjfdfdAdjfd - 6*
      tracefdAdjfdfdmH2I2Adjfd - 4*mHd2*tracefdAdjfdfuAdjfu - 4*mHu2*
      tracefdAdjfdfuAdjfu - 4*tracefdAdjfdfuAdjfuconjmSI2 - 4*
      tracefdAdjfdfumH1I2Adjfu - 12*tracefdAdjfdTfdAdjTfd - 4*
      tracefdAdjfdTfuAdjTfu - 12*tracefdAdjTfdTfdAdjfd - 4*
      tracefdAdjTfdTfuAdjfu - 2*tracefdconjLambda12TpLambda12AdjfdconjmSI2 - 2*
      tracefdconjTLambda12TpTLambda12Adjfd - 6*tracefdmH2I2AdjfdfdAdjfd - 4*
      tracefdmH2I2AdjfdfuAdjfu - 4*tracefuAdjfuTfdAdjTfd - 4*
      tracefuAdjTfuTfdAdjfd - 6*tracegDAdjgDconjmq2TpYdconjYd - 6*
      tracegDAdjgDTpTYdconjTYd - 6*tracegDAdjgDTpYdconjmd2conjYd - 6*mHd2*
      tracegDAdjgDTpYdconjYd - 6*mHp2*tracegDAdjgDTpYdconjYd - 6*
      tracegDAdjgDTpYdconjYdconjmq2 - 6*tracegDconjmDxbar2AdjgDTpYdconjYd - 4*
      tracehEAdjhEme2YeAdjYe - 4*tracehEAdjhETYeAdjTYe - 4*mHd2*
      tracehEAdjhEYeAdjYe - 4*mHp2*tracehEAdjhEYeAdjYe - 4*
      tracehEAdjhEYeAdjYeme2 - 4*tracehEAdjhEYeml2AdjYe - 4*
      tracehEAdjThETYeAdjYe - 4*tracehEmH1I2AdjhEYeAdjYe - 2*
      traceLambda12AdjLambda12conjmH2I2Tpfdconjfd - 2*mHd2*
      traceLambda12AdjLambda12Tpfdconjfd - 2*ms2*
      traceLambda12AdjLambda12Tpfdconjfd - 2*
      traceLambda12AdjLambda12TpfdconjfdconjmH2I2 - 2*
      traceLambda12AdjLambda12TpTfdconjTfd - 36*tracemd2YdAdjYdYdAdjYd - 6*
      tracemd2YdAdjYuYuAdjYd - 12*traceme2YeAdjYeYeAdjYe - 2*
      tracemH1I2AdjLambda12TpfdconjfdLambda12 - 12*traceml2AdjYeYeAdjYeYe - 36*
      tracemq2AdjYdYdAdjYdYd - 6*tracemq2AdjYdYdAdjYuYu - 6*
      tracemq2AdjYuYuAdjYdYd - 6*tracemu2YuAdjYdYdAdjYu - 36*
      traceYdAdjTYdTYdAdjYd - 6*traceYdAdjTYuTYuAdjYd - 36*
      traceYdAdjYdTYdAdjTYd - 36*mHd2*traceYdAdjYdYdAdjYd - 6*
      traceYdAdjYuTYuAdjTYd - 6*mHd2*traceYdAdjYuYuAdjYd - 6*mHu2*
      traceYdAdjYuYuAdjYd - 6*traceYdconjTgDTpTgDAdjYd - 4*
      traceYeAdjTYeThEAdjhE - 12*traceYeAdjTYeTYeAdjYe - 4*
      traceYeAdjYeThEAdjThE - 12*traceYeAdjYeTYeAdjTYe - 12*mHd2*
      traceYeAdjYeYeAdjYe - 6*traceYuAdjTYdTYdAdjYu - 6*traceYuAdjYdTYdAdjTYu +
      87*Power(g2,4)*AbsSqr(MassWB) - 2*traceconjTfuTpTfu*AbsSqr(Lambdax) - 6*
      traceconjTKappaTpTKappa*AbsSqr(Lambdax) - 4*traceconjTLambda12TpTLambda12
      *AbsSqr(Lambdax) - 6*traceconjTYuTpTYu*AbsSqr(Lambdax) - 2*mHd2*
      tracefuAdjfu*AbsSqr(Lambdax) - 4*mHu2*tracefuAdjfu*AbsSqr(Lambdax) - 2*
      ms2*tracefuAdjfu*AbsSqr(Lambdax) - 2*tracefuAdjfuconjmSI2*AbsSqr(Lambdax)
      - 2*tracefumH1I2Adjfu*AbsSqr(Lambdax) - 6*mHd2*traceKappaAdjKappa*AbsSqr
      (Lambdax) - 6*mHu2*traceKappaAdjKappa*AbsSqr(Lambdax) - 12*ms2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 6*traceKappaAdjKappaconjmDx2*AbsSqr(
      Lambdax) - 6*traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 4*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 4*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 8*ms2*traceLambda12AdjLambda12
      *AbsSqr(Lambdax) - 4*traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) -
      4*tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 6*tracemq2AdjYuYu*
      AbsSqr(Lambdax) - 6*tracemu2YuAdjYu*AbsSqr(Lambdax) - 6*mHd2*traceYuAdjYu
      *AbsSqr(Lambdax) - 12*mHu2*traceYuAdjYu*AbsSqr(Lambdax) - 6*ms2*
      traceYuAdjYu*AbsSqr(Lambdax) - 2*mHd2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*
      mHu2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*mphi2*AbsSqr(Lambdax)*AbsSqr(
      Sigmax) - 4*ms2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*msbar2*AbsSqr(Lambdax)
      *AbsSqr(Sigmax) - 2*tracefuAdjfu*AbsSqr(TLambdax) - 6*traceKappaAdjKappa*
      AbsSqr(TLambdax) - 4*traceLambda12AdjLambda12*AbsSqr(TLambdax) - 6*
      traceYuAdjYu*AbsSqr(TLambdax) - 24*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 2*
      AbsSqr(Sigmax)*AbsSqr(TLambdax) - 2*AbsSqr(Lambdax)*AbsSqr(TSigmax) - 2*
      traceAdjfuTfu*Conj(TLambdax)*Lambdax - 6*traceAdjKappaTKappa*Conj(
      TLambdax)*Lambdax - 4*traceAdjLambda12TLambda12*Conj(TLambdax)*Lambdax -
      6*traceAdjYuTYu*Conj(TLambdax)*Lambdax + 1.2*Tr2U111*Sqr(g1) - 0.8*
      traceconjTYdTpTYd*Sqr(g1) + 0.8*MassB*traceconjTYdTpYd*Sqr(g1) + 2.4*
      traceconjTYeTpTYe*Sqr(g1) - 2.4*MassB*traceconjTYeTpYe*Sqr(g1) - 0.8*
      tracemd2YdAdjYd*Sqr(g1) + 2.4*traceme2YeAdjYe*Sqr(g1) + 2.4*
      traceml2AdjYeYe*Sqr(g1) - 0.8*tracemq2AdjYdYd*Sqr(g1) - 0.8*mHd2*
      traceYdAdjYd*Sqr(g1) + 2.4*mHd2*traceYeAdjYe*Sqr(g1) + 1.8*Tr2U144*Sqr(
      g1p) - 2*MassBp*traceconjTfdTpfd*Sqr(g1p) + 2*traceconjTfdTpTfd*Sqr(g1p)
      - 1.2*traceconjTYdTpTYd*Sqr(g1p) + 1.2*MassBp*traceconjTYdTpYd*Sqr(g1p) -
      0.4*traceconjTYeTpTYe*Sqr(g1p) + 0.4*MassBp*traceconjTYeTpYe*Sqr(g1p) +
      2*mHd2*tracefdAdjfd*Sqr(g1p) + 2*tracefdAdjfdconjmSI2*Sqr(g1p) + 2*
      tracefdmH2I2Adjfd*Sqr(g1p) - 1.2*tracemd2YdAdjYd*Sqr(g1p) - 0.4*
      traceme2YeAdjYe*Sqr(g1p) - 0.4*traceml2AdjYeYe*Sqr(g1p) - 1.2*
      tracemq2AdjYdYd*Sqr(g1p) - 1.2*mHd2*traceYdAdjYd*Sqr(g1p) - 0.4*mHd2*
      traceYeAdjYe*Sqr(g1p) - 0.5*mHd2*AbsSqr(Lambdax)*Sqr(g1p) - 0.5*mHu2*
      AbsSqr(Lambdax)*Sqr(g1p) - 0.5*ms2*AbsSqr(Lambdax)*Sqr(g1p) - 0.5*AbsSqr(
      TLambdax)*Sqr(g1p) + 0.5*MassBp*Conj(TLambdax)*Lambdax*Sqr(g1p) + 3.6*
      AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) +
      5.4*AbsSqr(MassWB)*Sqr(g1p)*Sqr(g2) + 2.7*MassBp*Conj(MassWB)*Sqr(g1p)*
      Sqr(g2) + 0.02*Conj(MassB)*Sqr(g1)*(40*traceAdjYdTYd - 120*traceAdjYeTYe
      - 80*MassB*traceYdAdjYd + 240*MassB*traceYeAdjYe + 1782*MassB*Sqr(g1) -
      18*MassB*Sqr(g1p) - 9*MassBp*Sqr(g1p) + 180*MassB*Sqr(g2) + 90*MassWB*Sqr
      (g2)) + 32*traceconjTYdTpTYd*Sqr(g3) - 32*MassG*traceconjTYdTpYd*Sqr(g3)
      + 32*tracemd2YdAdjYd*Sqr(g3) + 32*tracemq2AdjYdYd*Sqr(g3) + 32*mHd2*
      traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) + 0.1*mHd2*AbsSqr(Lambdax)*Sqr(g1p)*Sqr
      (QS) + 0.1*mHu2*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS) + 0.1*ms2*AbsSqr(Lambdax
      )*Sqr(g1p)*Sqr(QS) + 0.1*AbsSqr(TLambdax)*Sqr(g1p)*Sqr(QS) - 0.1*MassBp*
      Conj(TLambdax)*Lambdax*Sqr(g1p)*Sqr(QS) - 12*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 12*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*ms2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 0.01*Conj(MassBp)*Sqr(g1p)*(-200*traceAdjfdTfd +
      120*traceAdjYdTYd + 40*traceAdjYeTYe + 400*MassBp*tracefdAdjfd - 240*
      MassBp*traceYdAdjYd - 80*MassBp*traceYeAdjYe - 18*MassB*Sqr(g1) - 36*
      MassBp*Sqr(g1) + 5319*MassBp*Sqr(g1p) + 540*MassBp*Sqr(g2) + 270*MassWB*
      Sqr(g2) + 27*MassBp*Sqr(g1p)*Sqr(QS) + 10*Conj(Lambdax)*(-5 + Sqr(QS))*(2
      *MassBp*Lambdax - TLambdax)) - 2*traceconjTfuTpfu*Conj(Lambdax)*TLambdax
      - 6*traceconjTKappaTpKappa*Conj(Lambdax)*TLambdax - 4*
      traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax - 6*traceconjTYuTpYu*
      Conj(Lambdax)*TLambdax - 2*Conj(Lambdax)*Conj(TSigmax)*Sigmax*TLambdax -
      2*Conj(Sigmax)*Conj(TLambdax)*Lambdax*TSigmax));


   return beta_mHd2;
}

} // namespace flexiblesusy
