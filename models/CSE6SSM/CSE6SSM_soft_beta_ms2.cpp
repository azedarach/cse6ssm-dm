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

// File generated at Wed 3 Jun 2015 23:43:31

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = Re(oneOver16PiSqr*(0.31622776601683794*g1p*QS*Tr14 + 6*
      traceconjTKappaTpTKappa + 4*traceconjTLambda12TpTLambda12 + 6*ms2*
      traceKappaAdjKappa + 6*traceKappaAdjKappaconjmDx2 + 6*
      traceKappaconjmDxbar2AdjKappa + 4*ms2*traceLambda12AdjLambda12 + 4*
      traceLambda12AdjLambda12conjmH2I2 + 4*tracemH1I2AdjLambda12Lambda12 + 4*(
      mHd2 + mHu2 + ms2)*AbsSqr(Lambdax) + 2*mphi2*AbsSqr(Sigmax) + 2*ms2*
      AbsSqr(Sigmax) + 2*msbar2*AbsSqr(Sigmax) + 4*AbsSqr(TLambdax) + 2*AbsSqr(
      TSigmax) - 0.2*AbsSqr(MassBp)*Sqr(g1p)*Sqr(QS)));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
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
   const double tracefdmH2I2Adjfd = TRACE_STRUCT.tracefdmH2I2Adjfd;
   const double tracefdAdjfdconjmSI2 = TRACE_STRUCT.tracefdAdjfdconjmSI2;
   const double tracefumH1I2Adjfu = TRACE_STRUCT.tracefumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2 = TRACE_STRUCT.tracefuAdjfuconjmSI2;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
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
   const double tracefdconjTLambda12TpTLambda12Adjfd =
      TRACE_STRUCT.tracefdconjTLambda12TpTLambda12Adjfd;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracefuAdjLambda12TLambda12AdjTfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12AdjTfu;
   const double tracefuAdjTLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjTLambda12TLambda12Adjfu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracegDAdjKappaTKappaAdjTgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjTgD;
   const double tracegDAdjTKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjTKappaTKappaAdjgD;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12TLambda12AdjThE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjThE;
   const double tracehEAdjTLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjTLambda12TLambda12AdjhE;
   const double traceKappaAdjgDTgDAdjTKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjTKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjTgDTgDAdjKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa;
   const double traceLambda12AdjfuTfuAdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjfuTfuAdjTLambda12;
   const double traceLambda12AdjhEThEAdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjTLambda12;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;
   const double traceLambda12AdjLambda12TpTfdconjTfd =
      TRACE_STRUCT.traceLambda12AdjLambda12TpTfdconjTfd;
   const double traceLambda12AdjTfuTfuAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTfuTfuAdjLambda12;
   const double traceLambda12AdjThEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjThEThEAdjLambda12;
   const double traceLambda12AdjTLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12;
   const double traceAdjfdTfdconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceAdjfdTfdconjTLambda12TpLambda12;
   const double traceAdjLambda12TpfdconjTfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjTfdTLambda12;
   const double tracefdconjLambda12TpLambda12AdjfdconjmSI2 =
      TRACE_STRUCT.tracefdconjLambda12TpLambda12AdjfdconjmSI2;
   const double tracefumH1I2AdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefumH1I2AdjLambda12Lambda12Adjfu;
   const double tracefuAdjLambda12Lambda12mH1I2Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12mH1I2Adjfu;
   const double tracefuAdjLambda12Lambda12AdjfuconjmSI2 =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12AdjfuconjmSI2;
   const double tracefuAdjLambda12conjmH2I2Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12conjmH2I2Lambda12Adjfu;
   const double tracegDAdjKappaKappaAdjgDconjmq2 =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgDconjmq2;
   const double tracegDAdjKappaKappaconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaconjmDxbar2AdjgD;
   const double tracegDAdjKappaconjmDx2KappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaconjmDx2KappaAdjgD;
   const double tracegDconjmDxbar2AdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjKappaKappaAdjgD;
   const double tracehEmH1I2AdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEmH1I2AdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12Lambda12mH1I2AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12mH1I2AdjhE;
   const double tracehEAdjLambda12Lambda12AdjhEme2 =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhEme2;
   const double tracehEAdjLambda12conjmH2I2Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12conjmH2I2Lambda12AdjhE;
   const double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
   const double tracemH1I2AdjLambda12TpfdconjfdLambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12TpfdconjfdLambda12;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
   const double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12conjmH2I2Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Tpfdconjfd;
   const double traceLambda12AdjLambda12TpfdconjfdconjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12TpfdconjfdconjmH2I2;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = Re(twoLoop*(1.2649110640673518*g1p*QS*Tr34 - 4*
      traceAdjfdTfdconjTLambda12TpLambda12 - 4*
      traceAdjLambda12TpfdconjTfdTLambda12 - 4*
      tracefdconjLambda12TpLambda12AdjfdconjmSI2 - 4*
      tracefdconjTLambda12TpTLambda12Adjfd - 4*
      tracefuAdjLambda12conjmH2I2Lambda12Adjfu - 4*mHu2*
      tracefuAdjLambda12Lambda12Adjfu - 4*ms2*tracefuAdjLambda12Lambda12Adjfu -
      4*tracefuAdjLambda12Lambda12AdjfuconjmSI2 - 4*
      tracefuAdjLambda12Lambda12mH1I2Adjfu - 4*
      tracefuAdjLambda12TLambda12AdjTfu - 4*tracefuAdjTLambda12TLambda12Adjfu -
      4*tracefumH1I2AdjLambda12Lambda12Adjfu - 12*
      tracegDAdjKappaconjmDx2KappaAdjgD - 12*mHp2*tracegDAdjKappaKappaAdjgD -
      12*ms2*tracegDAdjKappaKappaAdjgD - 12*tracegDAdjKappaKappaAdjgDconjmq2 -
      12*tracegDAdjKappaKappaconjmDxbar2AdjgD - 12*tracegDAdjKappaTKappaAdjTgD
      - 12*tracegDAdjTKappaTKappaAdjgD - 12*
      tracegDconjmDxbar2AdjKappaKappaAdjgD - 4*
      tracehEAdjLambda12conjmH2I2Lambda12AdjhE - 4*mHp2*
      tracehEAdjLambda12Lambda12AdjhE - 4*ms2*tracehEAdjLambda12Lambda12AdjhE -
      4*tracehEAdjLambda12Lambda12AdjhEme2 - 4*
      tracehEAdjLambda12Lambda12mH1I2AdjhE - 4*
      tracehEAdjLambda12TLambda12AdjThE - 4*tracehEAdjTLambda12TLambda12AdjhE -
      4*tracehEmH1I2AdjLambda12Lambda12AdjhE - 12*traceKappaAdjgDTgDAdjTKappa
      - 12*traceKappaAdjKappaconjmDx2KappaAdjKappa - 24*ms2*
      traceKappaAdjKappaKappaAdjKappa - 12*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 - 12*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa - 24*
      traceKappaAdjKappaTKappaAdjTKappa - 12*traceKappaAdjTgDTgDAdjKappa - 24*
      traceKappaAdjTKappaTKappaAdjKappa - 12*
      traceKappaconjmDxbar2AdjKappaKappaAdjKappa - 4*
      traceLambda12AdjfuTfuAdjTLambda12 - 4*traceLambda12AdjhEThEAdjTLambda12 -
      8*traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 - 4*
      traceLambda12AdjLambda12conjmH2I2Tpfdconjfd - 16*ms2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 8*
      traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 - 16*
      traceLambda12AdjLambda12TLambda12AdjTLambda12 - 4*mHd2*
      traceLambda12AdjLambda12Tpfdconjfd - 4*ms2*
      traceLambda12AdjLambda12Tpfdconjfd - 4*
      traceLambda12AdjLambda12TpfdconjfdconjmH2I2 - 4*
      traceLambda12AdjLambda12TpTfdconjTfd - 4*
      traceLambda12AdjTfuTfuAdjLambda12 - 4*traceLambda12AdjThEThEAdjLambda12 -
      16*traceLambda12AdjTLambda12TLambda12AdjLambda12 - 16*
      tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 - 4*
      tracemH1I2AdjLambda12TpfdconjfdLambda12 - 16*mphi2*AbsSqr(KappaPr)*AbsSqr
      (Sigmax) - 4*ms2*AbsSqr(KappaPr)*AbsSqr(Sigmax) - 4*msbar2*AbsSqr(KappaPr
      )*AbsSqr(Sigmax) - 4*mHp2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 4*mHpbar2*
      AbsSqr(Sigmax)*AbsSqr(SigmaL) - 8*mphi2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 4
      *ms2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 4*msbar2*AbsSqr(Sigmax)*AbsSqr(
      SigmaL) - 4*AbsSqr(Sigmax)*AbsSqr(TKappaPr) - 4*tracefdAdjfd*AbsSqr(
      TLambdax) - 4*tracefuAdjfu*AbsSqr(TLambdax) - 12*traceYdAdjYd*AbsSqr(
      TLambdax) - 4*traceYeAdjYe*AbsSqr(TLambdax) - 12*traceYuAdjYu*AbsSqr(
      TLambdax) - 4*AbsSqr(KappaPr)*AbsSqr(TSigmax) - 16*AbsSqr(Sigmax)*AbsSqr(
      TSigmax) - 4*AbsSqr(SigmaL)*AbsSqr(TSigmax) - 4*AbsSqr(Sigmax)*AbsSqr(
      TSigmaL) - 4*traceAdjfdTfd*Conj(TLambdax)*Lambdax - 4*traceAdjfuTfu*Conj(
      TLambdax)*Lambdax - 12*traceAdjYdTYd*Conj(TLambdax)*Lambdax - 4*
      traceAdjYeTYe*Conj(TLambdax)*Lambdax - 12*traceAdjYuTYu*Conj(TLambdax)*
      Lambdax - 1.6*MassB*traceconjTKappaTpKappa*Sqr(g1) + 1.6*
      traceconjTKappaTpTKappa*Sqr(g1) - 2.4*MassB*traceconjTLambda12TpLambda12*
      Sqr(g1) + 2.4*traceconjTLambda12TpTLambda12*Sqr(g1) + 1.6*ms2*
      traceKappaAdjKappa*Sqr(g1) + 1.6*traceKappaAdjKappaconjmDx2*Sqr(g1) + 1.6
      *traceKappaconjmDxbar2AdjKappa*Sqr(g1) + 2.4*ms2*traceLambda12AdjLambda12
      *Sqr(g1) + 2.4*traceLambda12AdjLambda12conjmH2I2*Sqr(g1) + 2.4*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1) + 3.2*traceKappaAdjKappa*AbsSqr(
      MassB)*Sqr(g1) + 4.8*traceLambda12AdjLambda12*AbsSqr(MassB)*Sqr(g1) + 2.4
      *AbsSqr(TLambdax)*Sqr(g1) - 1.6*traceAdjKappaTKappa*Conj(MassB)*Sqr(g1) -
      2.4*traceAdjLambda12TLambda12*Conj(MassB)*Sqr(g1) - 2.4*MassB*Conj(
      TLambdax)*Lambdax*Sqr(g1) - 3.9*MassBp*traceconjTKappaTpKappa*Sqr(g1p) +
      3.9*traceconjTKappaTpTKappa*Sqr(g1p) - 2.6*MassBp*
      traceconjTLambda12TpLambda12*Sqr(g1p) + 2.6*traceconjTLambda12TpTLambda12
      *Sqr(g1p) + 3.9*ms2*traceKappaAdjKappa*Sqr(g1p) + 3.9*
      traceKappaAdjKappaconjmDx2*Sqr(g1p) + 3.9*traceKappaconjmDxbar2AdjKappa*
      Sqr(g1p) + 2.6*ms2*traceLambda12AdjLambda12*Sqr(g1p) + 2.6*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g1p) + 2.6*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1p) + 2.6*AbsSqr(TLambdax)*Sqr(g1p) -
      2.6*MassBp*Conj(TLambdax)*Lambdax*Sqr(g1p) - 12*MassWB*
      traceconjTLambda12TpLambda12*Sqr(g2) + 12*traceconjTLambda12TpTLambda12*
      Sqr(g2) + 12*ms2*traceLambda12AdjLambda12*Sqr(g2) + 12*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g2) + 12*
      tracemH1I2AdjLambda12Lambda12*Sqr(g2) + 24*traceLambda12AdjLambda12*
      AbsSqr(MassWB)*Sqr(g2) + 12*AbsSqr(TLambdax)*Sqr(g2) - 12*
      traceAdjLambda12TLambda12*Conj(MassWB)*Sqr(g2) - 12*MassWB*Conj(TLambdax)
      *Lambdax*Sqr(g2) - 32*MassG*traceconjTKappaTpKappa*Sqr(g3) + 32*
      traceconjTKappaTpTKappa*Sqr(g3) + 32*ms2*traceKappaAdjKappa*Sqr(g3) + 32*
      traceKappaAdjKappaconjmDx2*Sqr(g3) + 32*traceKappaconjmDxbar2AdjKappa*Sqr
      (g3) + 64*traceKappaAdjKappa*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjKappaTKappa*Conj(MassG)*Sqr(g3) + 0.2*Tr2U144*Sqr(g1p)*Sqr(QS) +
      0.3*MassBp*traceconjTKappaTpKappa*Sqr(g1p)*Sqr(QS) - 0.3*
      traceconjTKappaTpTKappa*Sqr(g1p)*Sqr(QS) + 0.2*MassBp*
      traceconjTLambda12TpLambda12*Sqr(g1p)*Sqr(QS) - 0.2*
      traceconjTLambda12TpTLambda12*Sqr(g1p)*Sqr(QS) - 0.3*ms2*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) - 0.3*traceKappaAdjKappaconjmDx2*Sqr(
      g1p)*Sqr(QS) - 0.3*traceKappaconjmDxbar2AdjKappa*Sqr(g1p)*Sqr(QS) - 0.2*
      ms2*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS) - 0.2*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g1p)*Sqr(QS) - 0.2*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1p)*Sqr(QS) - 0.2*AbsSqr(TLambdax)*Sqr
      (g1p)*Sqr(QS) + 0.2*MassBp*Conj(TLambdax)*Lambdax*Sqr(g1p)*Sqr(QS) - 16*(
      mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 8*mphi2*Sqr(Conj(
      Sigmax))*Sqr(Sigmax) - 8*ms2*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 8*msbar2*Sqr
      (Conj(Sigmax))*Sqr(Sigmax) - 4*Conj(KappaPr)*Conj(TSigmax)*Sigmax*
      TKappaPr + 0.02*Conj(MassBp)*Sqr(g1p)*(-195*traceAdjKappaTKappa - 130*
      traceAdjLambda12TLambda12 + 260*MassBp*traceLambda12AdjLambda12 + 3*
      MassBp*Power(QS,4)*Sqr(g1p) - 30*MassBp*traceKappaAdjKappa*(-13 + Sqr(QS)
      ) + 15*traceAdjKappaTKappa*Sqr(QS) + 10*traceAdjLambda12TLambda12*Sqr(QS)
      - 20*MassBp*traceLambda12AdjLambda12*Sqr(QS) + 282*MassBp*Sqr(g1p)*Sqr(
      QS) - 10*Conj(Lambdax)*(-13 + Sqr(QS))*(2*MassBp*Lambdax - TLambdax)) +
      Conj(Lambdax)*(-4*traceconjTfdTpTfd*Lambdax - 4*traceconjTfuTpTfu*Lambdax
      - 12*traceconjTYdTpTYd*Lambdax - 4*traceconjTYeTpTYe*Lambdax - 12*
      traceconjTYuTpTYu*Lambdax - 8*mHd2*tracefdAdjfd*Lambdax - 4*mHu2*
      tracefdAdjfd*Lambdax - 4*ms2*tracefdAdjfd*Lambdax - 4*
      tracefdAdjfdconjmSI2*Lambdax - 4*tracefdmH2I2Adjfd*Lambdax - 4*mHd2*
      tracefuAdjfu*Lambdax - 8*mHu2*tracefuAdjfu*Lambdax - 4*ms2*tracefuAdjfu*
      Lambdax - 4*tracefuAdjfuconjmSI2*Lambdax - 4*tracefumH1I2Adjfu*Lambdax -
      12*tracemd2YdAdjYd*Lambdax - 4*traceme2YeAdjYe*Lambdax - 4*
      traceml2AdjYeYe*Lambdax - 12*tracemq2AdjYdYd*Lambdax - 12*tracemq2AdjYuYu
      *Lambdax - 12*tracemu2YuAdjYu*Lambdax - 24*mHd2*traceYdAdjYd*Lambdax - 12
      *mHu2*traceYdAdjYd*Lambdax - 12*ms2*traceYdAdjYd*Lambdax - 8*mHd2*
      traceYeAdjYe*Lambdax - 4*mHu2*traceYeAdjYe*Lambdax - 4*ms2*traceYeAdjYe*
      Lambdax - 12*mHd2*traceYuAdjYu*Lambdax - 24*mHu2*traceYuAdjYu*Lambdax -
      12*ms2*traceYuAdjYu*Lambdax - 32*AbsSqr(TLambdax)*Lambdax + 2.4*mHd2*
      Lambdax*Sqr(g1) + 2.4*mHu2*Lambdax*Sqr(g1) + 2.4*ms2*Lambdax*Sqr(g1) +
      2.6*mHd2*Lambdax*Sqr(g1p) + 2.6*mHu2*Lambdax*Sqr(g1p) + 2.6*ms2*Lambdax*
      Sqr(g1p) + 12*mHd2*Lambdax*Sqr(g2) + 12*mHu2*Lambdax*Sqr(g2) + 12*ms2*
      Lambdax*Sqr(g2) - 0.2*mHd2*Lambdax*Sqr(g1p)*Sqr(QS) - 0.2*mHu2*Lambdax*
      Sqr(g1p)*Sqr(QS) - 0.2*ms2*Lambdax*Sqr(g1p)*Sqr(QS) + 2.4*Conj(MassB)*Sqr
      (g1)*(2*MassB*Lambdax - TLambdax) + 12*Conj(MassWB)*Sqr(g2)*(2*MassWB*
      Lambdax - TLambdax) - 4*traceconjTfdTpfd*TLambdax - 4*traceconjTfuTpfu*
      TLambdax - 12*traceconjTYdTpYd*TLambdax - 4*traceconjTYeTpYe*TLambdax -
      12*traceconjTYuTpYu*TLambdax) - 4*Conj(Sigmax)*Conj(TKappaPr)*KappaPr*
      TSigmax - 4*Conj(Sigmax)*Conj(TSigmaL)*SigmaL*TSigmax - 4*Conj(SigmaL)*
      Conj(TSigmax)*Sigmax*TSigmaL));


   return beta_ms2;
}

} // namespace flexiblesusy
