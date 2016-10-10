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

// File generated at Wed 3 Jun 2015 23:43:42

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHpbar2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHpbar2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHpbar2;

   beta_mHpbar2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 -
      0.6324555320336759*g1p*Tr14 + 2*mHp2*AbsSqr(SigmaL) + 2*mHpbar2*AbsSqr(
      SigmaL) + 2*mphi2*AbsSqr(SigmaL) + 2*AbsSqr(TSigmaL) - 1.2*AbsSqr(MassB)*
      Sqr(g1) - 0.8*AbsSqr(MassBp)*Sqr(g1p) - 6*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHpbar2;
}

/**
 * Calculates the two-loop beta function of mHpbar2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_mHpbar2_two_loop(const Soft_traces& soft_traces) const
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
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHpbar2;

   beta_mHpbar2 = Re(twoLoop*(6*Power(g2,4)*Tr22 - 0.9797958971132712*g1*
      g1p*Tr2U114 - 0.9797958971132712*g1*g1p*Tr2U141 + 3.0983866769659336*g1*
      Tr31 - 2.5298221281347035*g1p*Tr34 + 87*Power(g2,4)*AbsSqr(MassWB) - 6*
      traceconjTgDTpTgD*AbsSqr(SigmaL) - 2*traceconjThETpThE*AbsSqr(SigmaL) -
      12*mHp2*tracegDAdjgD*AbsSqr(SigmaL) - 6*mHpbar2*tracegDAdjgD*AbsSqr(
      SigmaL) - 6*mphi2*tracegDAdjgD*AbsSqr(SigmaL) - 6*tracegDAdjgDconjmq2*
      AbsSqr(SigmaL) - 6*tracegDconjmDxbar2AdjgD*AbsSqr(SigmaL) - 4*mHp2*
      tracehEAdjhE*AbsSqr(SigmaL) - 2*mHpbar2*tracehEAdjhE*AbsSqr(SigmaL) - 2*
      mphi2*tracehEAdjhE*AbsSqr(SigmaL) - 2*tracehEAdjhEme2*AbsSqr(SigmaL) - 2*
      tracehEmH1I2AdjhE*AbsSqr(SigmaL) - 4*mHp2*AbsSqr(KappaPr)*AbsSqr(SigmaL)
      - 4*mHpbar2*AbsSqr(KappaPr)*AbsSqr(SigmaL) - 16*mphi2*AbsSqr(KappaPr)*
      AbsSqr(SigmaL) - 2*mHp2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 2*mHpbar2*AbsSqr(
      Sigmax)*AbsSqr(SigmaL) - 4*mphi2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 2*ms2*
      AbsSqr(Sigmax)*AbsSqr(SigmaL) - 2*msbar2*AbsSqr(Sigmax)*AbsSqr(SigmaL) -
      4*AbsSqr(SigmaL)*AbsSqr(TKappaPr) - 2*AbsSqr(SigmaL)*AbsSqr(TSigmax) - 6*
      tracegDAdjgD*AbsSqr(TSigmaL) - 2*tracehEAdjhE*AbsSqr(TSigmaL) - 4*AbsSqr(
      KappaPr)*AbsSqr(TSigmaL) - 2*AbsSqr(Sigmax)*AbsSqr(TSigmaL) - 24*AbsSqr(
      SigmaL)*AbsSqr(TSigmaL) - 6*traceAdjgDTgD*Conj(TSigmaL)*SigmaL - 2*
      traceAdjhEThE*Conj(TSigmaL)*SigmaL + 1.2*Tr2U111*Sqr(g1) + 0.8*Tr2U144*
      Sqr(g1p) + 3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 1.8*MassB*Conj(MassWB)*
      Sqr(g1)*Sqr(g2) + 2.4*AbsSqr(MassWB)*Sqr(g1p)*Sqr(g2) + 1.2*MassBp*Conj(
      MassWB)*Sqr(g1p)*Sqr(g2) + 0.36*Conj(MassB)*Sqr(g1)*(99*MassB*Sqr(g1) + 2
      *(2*MassB + MassBp)*Sqr(g1p) + 5*(2*MassB + MassWB)*Sqr(g2)) + 0.12*Conj(
      MassBp)*Sqr(g1p)*(6*(MassB + 2*MassBp)*Sqr(g1) + 10*(2*MassBp + MassWB)*
      Sqr(g2) + MassBp*Sqr(g1p)*(192 + Sqr(QS))) - 12*mHp2*Sqr(Conj(SigmaL))*
      Sqr(SigmaL) - 12*mHpbar2*Sqr(Conj(SigmaL))*Sqr(SigmaL) - 12*mphi2*Sqr(
      Conj(SigmaL))*Sqr(SigmaL) - 4*Conj(KappaPr)*Conj(TSigmaL)*SigmaL*TKappaPr
      - 2*Conj(Sigmax)*Conj(TSigmaL)*SigmaL*TSigmax - 6*traceconjTgDTpgD*Conj(
      SigmaL)*TSigmaL - 2*traceconjThETphE*Conj(SigmaL)*TSigmaL - 4*Conj(SigmaL
      )*Conj(TKappaPr)*KappaPr*TSigmaL - 2*Conj(SigmaL)*Conj(TSigmax)*Sigmax*
      TSigmaL));


   return beta_mHpbar2;
}

} // namespace flexiblesusy
