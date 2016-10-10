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

// File generated at Wed 3 Jun 2015 23:43:43

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mphi2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mphi2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_mphi2;

   beta_mphi2 = Re(2*oneOver16PiSqr*(6*mphi2*AbsSqr(KappaPr) + (mphi2 +
      ms2 + msbar2)*AbsSqr(Sigmax) + 2*mHp2*AbsSqr(SigmaL) + 2*mHpbar2*AbsSqr(
      SigmaL) + 2*mphi2*AbsSqr(SigmaL) + 2*AbsSqr(TKappaPr) + AbsSqr(TSigmax) +
      2*AbsSqr(TSigmaL)));


   return beta_mphi2;
}

/**
 * Calculates the two-loop beta function of mphi2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mphi2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDAdjgDconjmq2 = TRACE_STRUCT.tracegDAdjgDconjmq2;
   const double tracegDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgD;
   const double tracehEmH1I2AdjhE = TRACE_STRUCT.tracehEmH1I2AdjhE;
   const double tracehEAdjhEme2 = TRACE_STRUCT.tracehEAdjhEme2;
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
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;


   double beta_mphi2;

   beta_mphi2 = Re(0.2*twoLoop*(-480*mphi2*Sqr(Conj(KappaPr))*Sqr(KappaPr
      ) - 40*(mphi2 + ms2 + msbar2)*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 80*(mHp2 +
      mHpbar2 + mphi2)*Sqr(Conj(SigmaL))*Sqr(SigmaL) - 40*Conj(KappaPr)*((4*
      mphi2 + ms2 + msbar2)*AbsSqr(Sigmax)*KappaPr + 2*(mHp2 + mHpbar2 + 4*
      mphi2)*AbsSqr(SigmaL)*KappaPr + 8*AbsSqr(TKappaPr)*KappaPr + AbsSqr(
      TSigmax)*KappaPr + 2*AbsSqr(TSigmaL)*KappaPr + Conj(TSigmax)*Sigmax*
      TKappaPr + 2*Conj(TSigmaL)*SigmaL*TKappaPr) + Conj(Sigmax)*(-30*
      traceconjTKappaTpTKappa*Sigmax - 20*traceconjTLambda12TpTLambda12*Sigmax
      - 30*mphi2*traceKappaAdjKappa*Sigmax - 60*ms2*traceKappaAdjKappa*Sigmax -
      30*msbar2*traceKappaAdjKappa*Sigmax - 30*traceKappaAdjKappaconjmDx2*
      Sigmax - 30*traceKappaconjmDxbar2AdjKappa*Sigmax - 20*mphi2*
      traceLambda12AdjLambda12*Sigmax - 40*ms2*traceLambda12AdjLambda12*Sigmax
      - 20*msbar2*traceLambda12AdjLambda12*Sigmax - 20*
      traceLambda12AdjLambda12conjmH2I2*Sigmax - 20*
      tracemH1I2AdjLambda12Lambda12*Sigmax - 20*(mHd2 + mHu2 + mphi2 + 2*ms2 +
      msbar2)*AbsSqr(Lambdax)*Sigmax - 40*AbsSqr(TKappaPr)*Sigmax - 20*AbsSqr(
      TLambdax)*Sigmax - 80*AbsSqr(TSigmax)*Sigmax + mphi2*Sigmax*Sqr(g1p)*Sqr(
      QS) + ms2*Sigmax*Sqr(g1p)*Sqr(QS) + msbar2*Sigmax*Sqr(g1p)*Sqr(QS) + Conj
      (MassBp)*Sqr(g1p)*Sqr(QS)*(2*MassBp*Sigmax - TSigmax) - 30*
      traceconjTKappaTpKappa*TSigmax - 20*traceconjTLambda12TpLambda12*TSigmax
      - 40*Conj(TKappaPr)*KappaPr*TSigmax - 20*Conj(TLambdax)*Lambdax*TSigmax)
      - Conj(TSigmax)*(Sigmax*(30*traceAdjKappaTKappa + 20*
      traceAdjLambda12TLambda12 + MassBp*Sqr(g1p)*Sqr(QS)) + (30*
      traceKappaAdjKappa + 20*traceLambda12AdjLambda12 - Sqr(g1p)*Sqr(QS))*
      TSigmax + 20*Conj(Lambdax)*(Sigmax*TLambdax + Lambdax*TSigmax)) - 4*Conj(
      TSigmaL)*(SigmaL*(15*traceAdjgDTgD + 5*traceAdjhEThE + 3*MassB*Sqr(g1) +
      2*MassBp*Sqr(g1p) + 15*MassWB*Sqr(g2)) + (15*tracegDAdjgD + 5*
      tracehEAdjhE - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2))*TSigmaL) + 4*Conj(
      SigmaL)*(-15*traceconjTgDTpTgD*SigmaL - 5*traceconjThETpThE*SigmaL - 30*
      mHp2*tracegDAdjgD*SigmaL - 15*mHpbar2*tracegDAdjgD*SigmaL - 15*mphi2*
      tracegDAdjgD*SigmaL - 15*tracegDAdjgDconjmq2*SigmaL - 15*
      tracegDconjmDxbar2AdjgD*SigmaL - 10*mHp2*tracehEAdjhE*SigmaL - 5*mHpbar2*
      tracehEAdjhE*SigmaL - 5*mphi2*tracehEAdjhE*SigmaL - 5*tracehEAdjhEme2*
      SigmaL - 5*tracehEmH1I2AdjhE*SigmaL - 20*AbsSqr(TKappaPr)*SigmaL - 40*
      AbsSqr(TSigmaL)*SigmaL + 3*mHp2*SigmaL*Sqr(g1) + 3*mHpbar2*SigmaL*Sqr(g1)
      + 3*mphi2*SigmaL*Sqr(g1) + 2*mHp2*SigmaL*Sqr(g1p) + 2*mHpbar2*SigmaL*Sqr
      (g1p) + 2*mphi2*SigmaL*Sqr(g1p) + 15*mHp2*SigmaL*Sqr(g2) + 15*mHpbar2*
      SigmaL*Sqr(g2) + 15*mphi2*SigmaL*Sqr(g2) + 30*AbsSqr(MassWB)*SigmaL*Sqr(
      g2) + 3*Conj(MassB)*Sqr(g1)*(2*MassB*SigmaL - TSigmaL) + 2*Conj(MassBp)*
      Sqr(g1p)*(2*MassBp*SigmaL - TSigmaL) - 15*traceconjTgDTpgD*TSigmaL - 5*
      traceconjThETphE*TSigmaL - 20*Conj(TKappaPr)*KappaPr*TSigmaL - 15*Conj(
      MassWB)*Sqr(g2)*TSigmaL)));


   return beta_mphi2;
}

} // namespace flexiblesusy
