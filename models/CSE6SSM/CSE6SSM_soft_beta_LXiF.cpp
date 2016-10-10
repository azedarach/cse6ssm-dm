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

// File generated at Wed 3 Jun 2015 23:43:17

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of LXiF.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_LXiF_one_loop(const Soft_traces& soft_traces) const
{


   double beta_LXiF;

   beta_LXiF = Re(oneOver16PiSqr*(2*MuPhi*BMuPhi*Conj(KappaPr) - 4*MuPhi*
      BMuPr*Conj(SigmaL) + 4*mphi2*Conj(MuPhi)*KappaPr - 4*mHp2*Conj(MuPr)*
      SigmaL - 4*mHpbar2*Conj(MuPr)*SigmaL + 2*AbsSqr(KappaPr)*LXiF + AbsSqr(
      Sigmax)*LXiF + 2*AbsSqr(SigmaL)*LXiF + 2*Conj(BMuPhi)*TKappaPr + 4*Conj(
      KappaPr)*XiF*TKappaPr + 2*Conj(Sigmax)*XiF*TSigmax - 4*Conj(BMuPr)*
      TSigmaL + 4*Conj(SigmaL)*XiF*TSigmaL));


   return beta_LXiF;
}

/**
 * Calculates the two-loop beta function of LXiF.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_LXiF_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDmDxbar2AdjgD = TRACE_STRUCT.tracegDmDxbar2AdjgD;
   const double tracegDAdjgDmq2 = TRACE_STRUCT.tracegDAdjgDmq2;
   const double tracehEAdjhEconjme2 = TRACE_STRUCT.tracehEAdjhEconjme2;
   const double tracehEconjmH1I2AdjhE =
      TRACE_STRUCT.tracehEconjmH1I2AdjhE;


   double beta_LXiF;

   beta_LXiF = Re(twoLoop*(-4*MuPhi*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax) +
      2*AbsSqr(SigmaL))*BMuPhi*Conj(KappaPr) + 12*MuPhi*tracegDAdjgD*BMuPr*
      Conj(SigmaL) + 4*MuPhi*tracehEAdjhE*BMuPr*Conj(SigmaL) + 12*MuPhi*
      traceAdjgDTgD*Conj(SigmaL)*MuPr + 4*MuPhi*traceAdjhEThE*Conj(SigmaL)*MuPr
      - 6*traceAdjKappaTKappa*AbsSqr(Sigmax)*XiF - 4*traceAdjLambda12TLambda12
      *AbsSqr(Sigmax)*XiF - 12*traceAdjgDTgD*AbsSqr(SigmaL)*XiF - 4*
      traceAdjhEThE*AbsSqr(SigmaL)*XiF + 12*traceAdjgDTgD*Conj(BMuPr)*SigmaL +
      4*traceAdjhEThE*Conj(BMuPr)*SigmaL + 12*traceconjTgDTpTgD*Conj(MuPr)*
      SigmaL + 4*traceconjThETpThE*Conj(MuPr)*SigmaL + 24*mHp2*tracegDAdjgD*
      Conj(MuPr)*SigmaL + 12*mHpbar2*tracegDAdjgD*Conj(MuPr)*SigmaL + 12*
      tracegDAdjgDmq2*Conj(MuPr)*SigmaL + 12*tracegDmDxbar2AdjgD*Conj(MuPr)*
      SigmaL + 8*mHp2*tracehEAdjhE*Conj(MuPr)*SigmaL + 4*mHpbar2*tracehEAdjhE*
      Conj(MuPr)*SigmaL + 4*tracehEAdjhEconjme2*Conj(MuPr)*SigmaL + 4*
      tracehEconjmH1I2AdjhE*Conj(MuPr)*SigmaL + 16*AbsSqr(TSigmaL)*Conj(MuPr)*
      SigmaL - 3*traceKappaAdjKappa*AbsSqr(Sigmax)*LXiF - 2*
      traceLambda12AdjLambda12*AbsSqr(Sigmax)*LXiF - 4*AbsSqr(KappaPr)*AbsSqr(
      Sigmax)*LXiF - 2*AbsSqr(Lambdax)*AbsSqr(Sigmax)*LXiF - 6*tracegDAdjgD*
      AbsSqr(SigmaL)*LXiF - 2*tracehEAdjhE*AbsSqr(SigmaL)*LXiF - 8*AbsSqr(
      KappaPr)*AbsSqr(SigmaL)*LXiF - 2.4*MuPhi*BMuPr*Conj(SigmaL)*Sqr(g1) + 2.4
      *MassB*MuPhi*Conj(SigmaL)*MuPr*Sqr(g1) - 2.4*MassB*AbsSqr(SigmaL)*XiF*Sqr
      (g1) + 2.4*MassB*Conj(BMuPr)*SigmaL*Sqr(g1) - 2.4*mHp2*Conj(MuPr)*SigmaL*
      Sqr(g1) - 2.4*mHpbar2*Conj(MuPr)*SigmaL*Sqr(g1) - 4.8*AbsSqr(MassB)*Conj(
      MuPr)*SigmaL*Sqr(g1) + 1.2*AbsSqr(SigmaL)*LXiF*Sqr(g1) - 1.6*MuPhi*BMuPr*
      Conj(SigmaL)*Sqr(g1p) + 1.6*MassBp*MuPhi*Conj(SigmaL)*MuPr*Sqr(g1p) - 1.6
      *MassBp*AbsSqr(SigmaL)*XiF*Sqr(g1p) + 1.6*MassBp*Conj(BMuPr)*SigmaL*Sqr(
      g1p) - 1.6*mHp2*Conj(MuPr)*SigmaL*Sqr(g1p) - 1.6*mHpbar2*Conj(MuPr)*
      SigmaL*Sqr(g1p) - 3.2*AbsSqr(MassBp)*Conj(MuPr)*SigmaL*Sqr(g1p) + 0.8*
      AbsSqr(SigmaL)*LXiF*Sqr(g1p) - 12*MuPhi*BMuPr*Conj(SigmaL)*Sqr(g2) + 12*
      MassWB*MuPhi*Conj(SigmaL)*MuPr*Sqr(g2) - 12*MassWB*AbsSqr(SigmaL)*XiF*Sqr
      (g2) + 12*MassWB*Conj(BMuPr)*SigmaL*Sqr(g2) - 12*mHp2*Conj(MuPr)*SigmaL*
      Sqr(g2) - 12*mHpbar2*Conj(MuPr)*SigmaL*Sqr(g2) - 24*AbsSqr(MassWB)*Conj(
      MuPr)*SigmaL*Sqr(g2) + 6*AbsSqr(SigmaL)*LXiF*Sqr(g2) - 0.2*MassBp*AbsSqr(
      Sigmax)*XiF*Sqr(g1p)*Sqr(QS) + 0.1*AbsSqr(Sigmax)*LXiF*Sqr(g1p)*Sqr(QS) +
      8*MuPhi*BMuPr*SigmaL*Sqr(Conj(SigmaL)) - 8*LXiF*Sqr(Conj(KappaPr))*Sqr(
      KappaPr) - 2*LXiF*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 16*mHp2*Conj(MuPr)*Conj
      (SigmaL)*Sqr(SigmaL) + 16*mHpbar2*Conj(MuPr)*Conj(SigmaL)*Sqr(SigmaL) + 8
      *mphi2*Conj(MuPr)*Conj(SigmaL)*Sqr(SigmaL) - 4*LXiF*Sqr(Conj(SigmaL))*Sqr
      (SigmaL) - 16*AbsSqr(KappaPr)*Conj(BMuPhi)*TKappaPr - 4*AbsSqr(Sigmax)*
      Conj(BMuPhi)*TKappaPr - 8*AbsSqr(SigmaL)*Conj(BMuPhi)*TKappaPr - 8*AbsSqr
      (Sigmax)*Conj(KappaPr)*XiF*TKappaPr - 16*AbsSqr(SigmaL)*Conj(KappaPr)*XiF
      *TKappaPr - 32*KappaPr*XiF*Sqr(Conj(KappaPr))*TKappaPr - 8*Sqr(MuPhi)*Sqr
      (Conj(KappaPr))*TKappaPr - 4*Conj(MuPhi)*((3*mphi2 + ms2 + msbar2)*AbsSqr
      (Sigmax)*KappaPr + 2*mHp2*AbsSqr(SigmaL)*KappaPr + 2*mHpbar2*AbsSqr(
      SigmaL)*KappaPr + 6*mphi2*AbsSqr(SigmaL)*KappaPr + 4*AbsSqr(TKappaPr)*
      KappaPr + AbsSqr(TSigmax)*KappaPr + 2*AbsSqr(TSigmaL)*KappaPr + 10*mphi2*
      Conj(KappaPr)*Sqr(KappaPr) + Conj(TSigmax)*Sigmax*TKappaPr + 2*Conj(
      TSigmaL)*SigmaL*TKappaPr) - 4*AbsSqr(Sigmax)*Conj(Lambdax)*XiF*TLambdax -
      4*Conj(BMuPhi)*Conj(Sigmax)*KappaPr*TSigmax - 6*traceKappaAdjKappa*Conj(
      Sigmax)*XiF*TSigmax - 4*traceLambda12AdjLambda12*Conj(Sigmax)*XiF*TSigmax
      - 8*AbsSqr(KappaPr)*Conj(Sigmax)*XiF*TSigmax - 4*AbsSqr(Lambdax)*Conj(
      Sigmax)*XiF*TSigmax - 4*Conj(KappaPr)*Conj(Sigmax)*Sqr(MuPhi)*TSigmax +
      0.2*Conj(Sigmax)*XiF*Sqr(g1p)*Sqr(QS)*TSigmax - 8*XiF*Sigmax*Sqr(Conj(
      Sigmax))*TSigmax + 12*tracegDAdjgD*Conj(BMuPr)*TSigmaL + 4*tracehEAdjhE*
      Conj(BMuPr)*TSigmaL + 16*AbsSqr(SigmaL)*Conj(BMuPr)*TSigmaL + 12*
      traceconjTgDTpgD*Conj(MuPr)*TSigmaL + 4*traceconjThETphE*Conj(MuPr)*
      TSigmaL - 8*Conj(BMuPhi)*Conj(SigmaL)*KappaPr*TSigmaL - 12*tracegDAdjgD*
      Conj(SigmaL)*XiF*TSigmaL - 4*tracehEAdjhE*Conj(SigmaL)*XiF*TSigmaL - 16*
      AbsSqr(KappaPr)*Conj(SigmaL)*XiF*TSigmaL - 2.4*Conj(BMuPr)*Sqr(g1)*
      TSigmaL + 2.4*MassB*Conj(MuPr)*Sqr(g1)*TSigmaL + 2.4*Conj(SigmaL)*XiF*Sqr
      (g1)*TSigmaL - 1.6*Conj(BMuPr)*Sqr(g1p)*TSigmaL + 1.6*MassBp*Conj(MuPr)*
      Sqr(g1p)*TSigmaL + 1.6*Conj(SigmaL)*XiF*Sqr(g1p)*TSigmaL - 12*Conj(BMuPr)
      *Sqr(g2)*TSigmaL + 12*MassWB*Conj(MuPr)*Sqr(g2)*TSigmaL + 12*Conj(SigmaL)
      *XiF*Sqr(g2)*TSigmaL - 8*Conj(KappaPr)*Conj(SigmaL)*Sqr(MuPhi)*TSigmaL +
      8*MuPhi*MuPr*Sqr(Conj(SigmaL))*TSigmaL - 16*XiF*SigmaL*Sqr(Conj(SigmaL))*
      TSigmaL));


   return beta_LXiF;
}

} // namespace flexiblesusy
