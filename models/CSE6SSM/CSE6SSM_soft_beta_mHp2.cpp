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

// File generated at Wed 3 Jun 2015 23:43:41

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHp2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHp2_one_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDAdjgDconjmq2 = TRACE_STRUCT.tracegDAdjgDconjmq2;
   const double tracegDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgD;
   const double tracehEmH1I2AdjhE = TRACE_STRUCT.tracehEmH1I2AdjhE;
   const double tracehEAdjhEme2 = TRACE_STRUCT.tracehEAdjhEme2;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHp2;

   beta_mHp2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 +
      0.6324555320336759*g1p*Tr14 + 6*traceconjTgDTpTgD + 2*traceconjThETpThE +
      6*mHp2*tracegDAdjgD + 6*tracegDAdjgDconjmq2 + 6*tracegDconjmDxbar2AdjgD
      + 2*mHp2*tracehEAdjhE + 2*tracehEAdjhEme2 + 2*tracehEmH1I2AdjhE + 2*mHp2*
      AbsSqr(SigmaL) + 2*mHpbar2*AbsSqr(SigmaL) + 2*mphi2*AbsSqr(SigmaL) + 2*
      AbsSqr(TSigmaL) - 1.2*AbsSqr(MassB)*Sqr(g1) - 0.8*AbsSqr(MassBp)*Sqr(g1p)
      - 6*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHp2;
}

/**
 * Calculates the two-loop beta function of mHp2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHp2_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDAdjgDconjmq2 = TRACE_STRUCT.tracegDAdjgDconjmq2;
   const double tracegDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgD;
   const double tracehEmH1I2AdjhE = TRACE_STRUCT.tracehEmH1I2AdjhE;
   const double tracehEAdjhEme2 = TRACE_STRUCT.tracehEAdjhEme2;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjhEThEAdjTfu =
      TRACE_STRUCT.tracefuAdjhEThEAdjTfu;
   const double tracefuAdjThEThEAdjfu =
      TRACE_STRUCT.tracefuAdjThEThEAdjfu;
   const double tracegDAdjgDgDAdjgD = TRACE_STRUCT.tracegDAdjgDgDAdjgD;
   const double tracegDAdjgDTgDAdjTgD =
      TRACE_STRUCT.tracegDAdjgDTgDAdjTgD;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjgDTpTYdconjTYd =
      TRACE_STRUCT.tracegDAdjgDTpTYdconjTYd;
   const double tracegDAdjgDTpTYuconjTYu =
      TRACE_STRUCT.tracegDAdjgDTpTYuconjTYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracegDAdjKappaTKappaAdjTgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjTgD;
   const double tracegDAdjTgDTgDAdjgD =
      TRACE_STRUCT.tracegDAdjTgDTgDAdjgD;
   const double tracegDAdjTKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjTKappaTKappaAdjgD;
   const double tracehEAdjfuTfuAdjThE =
      TRACE_STRUCT.tracehEAdjfuTfuAdjThE;
   const double tracehEAdjhEhEAdjhE = TRACE_STRUCT.tracehEAdjhEhEAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjhEThEAdjThE =
      TRACE_STRUCT.tracehEAdjhEThEAdjThE;
   const double tracehEAdjhETYeAdjTYe =
      TRACE_STRUCT.tracehEAdjhETYeAdjTYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double tracehEAdjLambda12TLambda12AdjThE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjThE;
   const double tracehEAdjTfuTfuAdjhE =
      TRACE_STRUCT.tracehEAdjTfuTfuAdjhE;
   const double tracehEAdjThEThEAdjhE =
      TRACE_STRUCT.tracehEAdjThEThEAdjhE;
   const double tracehEAdjThETYeAdjYe =
      TRACE_STRUCT.tracehEAdjThETYeAdjYe;
   const double tracehEAdjTLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjTLambda12TLambda12AdjhE;
   const double traceYdconjTgDTpTgDAdjYd =
      TRACE_STRUCT.traceYdconjTgDTpTgDAdjYd;
   const double traceYeAdjYeThEAdjThE =
      TRACE_STRUCT.traceYeAdjYeThEAdjThE;
   const double traceYeAdjTYeThEAdjhE =
      TRACE_STRUCT.traceYeAdjTYeThEAdjhE;
   const double traceYuconjTgDTpTgDAdjYu =
      TRACE_STRUCT.traceYuconjTgDTpTgDAdjYu;
   const double traceKappaAdjgDTgDAdjTKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjTKappa;
   const double traceKappaAdjTgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjTgDTgDAdjKappa;
   const double traceLambda12AdjhEThEAdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjTLambda12;
   const double traceLambda12AdjThEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjThEThEAdjLambda12;
   const double traceAdjgDTpYdconjTYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjTYdTgD;
   const double traceAdjgDTpYuconjTYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjTYuTgD;
   const double traceAdjYdTYdconjTgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjTgDTpgD;
   const double traceAdjYuTYuconjTgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjTgDTpgD;
   const double tracefumH1I2AdjhEhEAdjfu =
      TRACE_STRUCT.tracefumH1I2AdjhEhEAdjfu;
   const double tracefuAdjhEhEmH1I2Adjfu =
      TRACE_STRUCT.tracefuAdjhEhEmH1I2Adjfu;
   const double tracefuAdjhEhEAdjfuconjmSI2 =
      TRACE_STRUCT.tracefuAdjhEhEAdjfuconjmSI2;
   const double tracefuAdjhEme2hEAdjfu =
      TRACE_STRUCT.tracefuAdjhEme2hEAdjfu;
   const double tracegDAdjgDgDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDAdjgDgDconjmDxbar2AdjgD;
   const double tracegDAdjgDconjmq2gDAdjgD =
      TRACE_STRUCT.tracegDAdjgDconjmq2gDAdjgD;
   const double tracegDAdjgDconjmq2TpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDconjmq2TpYdconjYd;
   const double tracegDAdjgDconjmq2TpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDconjmq2TpYuconjYu;
   const double tracegDAdjgDTpYdconjmd2conjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjmd2conjYd;
   const double tracegDAdjgDTpYdconjYdconjmq2 =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYdconjmq2;
   const double tracegDAdjgDTpYuconjmu2conjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjmu2conjYu;
   const double tracegDAdjgDTpYuconjYuconjmq2 =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYuconjmq2;
   const double tracegDAdjKappaKappaAdjgDconjmq2 =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgDconjmq2;
   const double tracegDAdjKappaKappaconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaconjmDxbar2AdjgD;
   const double tracegDAdjKappaconjmDx2KappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaconjmDx2KappaAdjgD;
   const double tracegDconjmDxbar2AdjgDgDAdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgDgDAdjgD;
   const double tracegDconjmDxbar2AdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYdconjYd;
   const double tracegDconjmDxbar2AdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYuconjYu;
   const double tracegDconjmDxbar2AdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjKappaKappaAdjgD;
   const double tracehEmH1I2AdjhEhEAdjhE =
      TRACE_STRUCT.tracehEmH1I2AdjhEhEAdjhE;
   const double tracehEmH1I2AdjhEYeAdjYe =
      TRACE_STRUCT.tracehEmH1I2AdjhEYeAdjYe;
   const double tracehEmH1I2AdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEmH1I2AdjLambda12Lambda12AdjhE;
   const double tracehEAdjhEhEmH1I2AdjhE =
      TRACE_STRUCT.tracehEAdjhEhEmH1I2AdjhE;
   const double tracehEAdjhEhEAdjhEme2 =
      TRACE_STRUCT.tracehEAdjhEhEAdjhEme2;
   const double tracehEAdjhEme2hEAdjhE =
      TRACE_STRUCT.tracehEAdjhEme2hEAdjhE;
   const double tracehEAdjhEme2YeAdjYe =
      TRACE_STRUCT.tracehEAdjhEme2YeAdjYe;
   const double tracehEAdjhEYeml2AdjYe =
      TRACE_STRUCT.tracehEAdjhEYeml2AdjYe;
   const double tracehEAdjhEYeAdjYeme2 =
      TRACE_STRUCT.tracehEAdjhEYeAdjYeme2;
   const double tracehEAdjLambda12Lambda12mH1I2AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12mH1I2AdjhE;
   const double tracehEAdjLambda12Lambda12AdjhEme2 =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhEme2;
   const double tracehEAdjLambda12conjmH2I2Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12conjmH2I2Lambda12AdjhE;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHp2;

   beta_mHp2 = Re(twoLoop*(6*Power(g2,4)*Tr22 - 0.9797958971132712*g1*g1p
      *Tr2U114 - 0.9797958971132712*g1*g1p*Tr2U141 - 3.0983866769659336*g1*Tr31
      + 2.5298221281347035*g1p*Tr34 - 6*traceAdjgDTpYdconjTYdTgD - 6*
      traceAdjgDTpYuconjTYuTgD - 6*traceAdjYdTYdconjTgDTpgD - 6*
      traceAdjYuTYuconjTgDTpgD - 2*mHp2*tracefuAdjhEhEAdjfu - 2*mHu2*
      tracefuAdjhEhEAdjfu - 2*tracefuAdjhEhEAdjfuconjmSI2 - 2*
      tracefuAdjhEhEmH1I2Adjfu - 2*tracefuAdjhEme2hEAdjfu - 2*
      tracefuAdjhEThEAdjTfu - 2*tracefuAdjThEThEAdjfu - 2*
      tracefumH1I2AdjhEhEAdjfu - 36*tracegDAdjgDconjmq2gDAdjgD - 6*
      tracegDAdjgDconjmq2TpYdconjYd - 6*tracegDAdjgDconjmq2TpYuconjYu - 36*mHp2
      *tracegDAdjgDgDAdjgD - 18*tracegDAdjgDgDconjmDxbar2AdjgD - 36*
      tracegDAdjgDTgDAdjTgD - 6*tracegDAdjgDTpTYdconjTYd - 6*
      tracegDAdjgDTpTYuconjTYu - 6*tracegDAdjgDTpYdconjmd2conjYd - 6*mHd2*
      tracegDAdjgDTpYdconjYd - 6*mHp2*tracegDAdjgDTpYdconjYd - 6*
      tracegDAdjgDTpYdconjYdconjmq2 - 6*tracegDAdjgDTpYuconjmu2conjYu - 6*mHp2*
      tracegDAdjgDTpYuconjYu - 6*mHu2*tracegDAdjgDTpYuconjYu - 6*
      tracegDAdjgDTpYuconjYuconjmq2 - 6*tracegDAdjKappaconjmDx2KappaAdjgD - 6*
      mHp2*tracegDAdjKappaKappaAdjgD - 6*ms2*tracegDAdjKappaKappaAdjgD - 6*
      tracegDAdjKappaKappaAdjgDconjmq2 - 6*tracegDAdjKappaKappaconjmDxbar2AdjgD
      - 6*tracegDAdjKappaTKappaAdjTgD - 36*tracegDAdjTgDTgDAdjgD - 6*
      tracegDAdjTKappaTKappaAdjgD - 18*tracegDconjmDxbar2AdjgDgDAdjgD - 6*
      tracegDconjmDxbar2AdjgDTpYdconjYd - 6*tracegDconjmDxbar2AdjgDTpYuconjYu -
      6*tracegDconjmDxbar2AdjKappaKappaAdjgD - 2*tracehEAdjfuTfuAdjThE - 12*
      mHp2*tracehEAdjhEhEAdjhE - 6*tracehEAdjhEhEAdjhEme2 - 6*
      tracehEAdjhEhEmH1I2AdjhE - 6*tracehEAdjhEme2hEAdjhE - 4*
      tracehEAdjhEme2YeAdjYe - 12*tracehEAdjhEThEAdjThE - 4*
      tracehEAdjhETYeAdjTYe - 4*mHd2*tracehEAdjhEYeAdjYe - 4*mHp2*
      tracehEAdjhEYeAdjYe - 4*tracehEAdjhEYeAdjYeme2 - 4*tracehEAdjhEYeml2AdjYe
      - 2*tracehEAdjLambda12conjmH2I2Lambda12AdjhE - 2*mHp2*
      tracehEAdjLambda12Lambda12AdjhE - 2*ms2*tracehEAdjLambda12Lambda12AdjhE -
      2*tracehEAdjLambda12Lambda12AdjhEme2 - 2*
      tracehEAdjLambda12Lambda12mH1I2AdjhE - 2*
      tracehEAdjLambda12TLambda12AdjThE - 2*tracehEAdjTfuTfuAdjhE - 12*
      tracehEAdjThEThEAdjhE - 4*tracehEAdjThETYeAdjYe - 2*
      tracehEAdjTLambda12TLambda12AdjhE - 6*tracehEmH1I2AdjhEhEAdjhE - 4*
      tracehEmH1I2AdjhEYeAdjYe - 2*tracehEmH1I2AdjLambda12Lambda12AdjhE - 6*
      traceKappaAdjgDTgDAdjTKappa - 6*traceKappaAdjTgDTgDAdjKappa - 2*
      traceLambda12AdjhEThEAdjTLambda12 - 2*traceLambda12AdjThEThEAdjLambda12 -
      6*traceYdconjTgDTpTgDAdjYd - 4*traceYeAdjTYeThEAdjhE - 4*
      traceYeAdjYeThEAdjThE - 6*traceYuconjTgDTpTgDAdjYu + 87*Power(g2,4)*
      AbsSqr(MassWB) - 4*mHp2*AbsSqr(KappaPr)*AbsSqr(SigmaL) - 4*mHpbar2*AbsSqr
      (KappaPr)*AbsSqr(SigmaL) - 16*mphi2*AbsSqr(KappaPr)*AbsSqr(SigmaL) - 2*
      mHp2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 2*mHpbar2*AbsSqr(Sigmax)*AbsSqr(
      SigmaL) - 4*mphi2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 2*ms2*AbsSqr(Sigmax)*
      AbsSqr(SigmaL) - 2*msbar2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 4*AbsSqr(SigmaL
      )*AbsSqr(TKappaPr) - 2*AbsSqr(SigmaL)*AbsSqr(TSigmax) - 4*AbsSqr(KappaPr)
      *AbsSqr(TSigmaL) - 2*AbsSqr(Sigmax)*AbsSqr(TSigmaL) - 24*AbsSqr(SigmaL)*
      AbsSqr(TSigmaL) + 1.2*Tr2U111*Sqr(g1) + 0.8*MassB*traceconjTgDTpgD*Sqr(g1
      ) - 0.8*traceconjTgDTpTgD*Sqr(g1) - 2.4*MassB*traceconjThETphE*Sqr(g1) +
      2.4*traceconjThETpThE*Sqr(g1) - 0.8*mHp2*tracegDAdjgD*Sqr(g1) - 0.8*
      tracegDAdjgDconjmq2*Sqr(g1) - 0.8*tracegDconjmDxbar2AdjgD*Sqr(g1) + 2.4*
      mHp2*tracehEAdjhE*Sqr(g1) + 2.4*tracehEAdjhEme2*Sqr(g1) + 2.4*
      tracehEmH1I2AdjhE*Sqr(g1) + 0.8*Tr2U144*Sqr(g1p) - 1.8*MassBp*
      traceconjTgDTpgD*Sqr(g1p) + 1.8*traceconjTgDTpTgD*Sqr(g1p) - 0.6*MassBp*
      traceconjThETphE*Sqr(g1p) + 0.6*traceconjThETpThE*Sqr(g1p) + 1.8*mHp2*
      tracegDAdjgD*Sqr(g1p) + 1.8*tracegDAdjgDconjmq2*Sqr(g1p) + 1.8*
      tracegDconjmDxbar2AdjgD*Sqr(g1p) + 0.6*mHp2*tracehEAdjhE*Sqr(g1p) + 0.6*
      tracehEAdjhEme2*Sqr(g1p) + 0.6*tracehEmH1I2AdjhE*Sqr(g1p) + 3.6*AbsSqr(
      MassWB)*Sqr(g1)*Sqr(g2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 2.4*
      AbsSqr(MassWB)*Sqr(g1p)*Sqr(g2) + 1.2*MassBp*Conj(MassWB)*Sqr(g1p)*Sqr(g2
      ) + 0.04*Conj(MassB)*Sqr(g1)*(20*traceAdjgDTgD - 60*traceAdjhEThE - 40*
      MassB*tracegDAdjgD + 120*MassB*tracehEAdjhE + 891*MassB*Sqr(g1) + 36*
      MassB*Sqr(g1p) + 18*MassBp*Sqr(g1p) + 90*MassB*Sqr(g2) + 45*MassWB*Sqr(g2
      )) - 32*MassG*traceconjTgDTpgD*Sqr(g3) + 32*traceconjTgDTpTgD*Sqr(g3) +
      32*mHp2*tracegDAdjgD*Sqr(g3) + 32*tracegDAdjgDconjmq2*Sqr(g3) + 32*
      tracegDconjmDxbar2AdjgD*Sqr(g3) + 64*tracegDAdjgD*AbsSqr(MassG)*Sqr(g3) -
      32*traceAdjgDTgD*Conj(MassG)*Sqr(g3) + 0.12*Conj(MassBp)*Sqr(g1p)*(-15*
      traceAdjgDTgD - 5*traceAdjhEThE + 30*MassBp*tracegDAdjgD + 10*MassBp*
      tracehEAdjhE + 6*MassB*Sqr(g1) + 12*MassBp*Sqr(g1) + 192*MassBp*Sqr(g1p)
      + 20*MassBp*Sqr(g2) + 10*MassWB*Sqr(g2) + MassBp*Sqr(g1p)*Sqr(QS)) - 12*
      mHp2*Sqr(Conj(SigmaL))*Sqr(SigmaL) - 12*mHpbar2*Sqr(Conj(SigmaL))*Sqr(
      SigmaL) - 12*mphi2*Sqr(Conj(SigmaL))*Sqr(SigmaL) - 4*Conj(KappaPr)*Conj(
      TSigmaL)*SigmaL*TKappaPr - 2*Conj(Sigmax)*Conj(TSigmaL)*SigmaL*TSigmax -
      4*Conj(SigmaL)*Conj(TKappaPr)*KappaPr*TSigmaL - 2*Conj(SigmaL)*Conj(
      TSigmax)*Sigmax*TSigmaL));


   return beta_mHp2;
}

} // namespace flexiblesusy
