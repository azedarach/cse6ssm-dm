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

// File generated at Wed 3 Jun 2015 23:43:15

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of BMuPr.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_BMuPr_one_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;


   double beta_BMuPr;

   beta_BMuPr = Re(oneOver16PiSqr*(-2*BMuPhi*Conj(KappaPr)*SigmaL + BMuPr
      *(3*tracegDAdjgD + tracehEAdjhE + 6*AbsSqr(SigmaL) - 0.6*Sqr(g1) - 0.4*
      Sqr(g1p) - 3*Sqr(g2)) + 0.4*MuPr*(15*traceAdjgDTgD + 5*traceAdjhEThE + 3*
      MassB*Sqr(g1) + 2*MassBp*Sqr(g1p) + 15*MassWB*Sqr(g2) + 10*Conj(SigmaL)*
      TSigmaL)));


   return beta_BMuPr;
}

/**
 * Calculates the two-loop beta function of BMuPr.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_BMuPr_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracefuAdjhEThEAdjfu = TRACE_STRUCT.tracefuAdjhEThEAdjfu;
   const double tracegDAdjgDTgDAdjgD = TRACE_STRUCT.tracegDAdjgDTgDAdjgD;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjfuTfuAdjhE = TRACE_STRUCT.tracehEAdjfuTfuAdjhE;
   const double tracehEAdjhEThEAdjhE = TRACE_STRUCT.tracehEAdjhEThEAdjhE;
   const double tracehEAdjhETYeAdjYe = TRACE_STRUCT.tracehEAdjhETYeAdjYe;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracegDAdjgDgDAdjgD = TRACE_STRUCT.tracegDAdjgDgDAdjgD;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjhEhEAdjhE = TRACE_STRUCT.tracehEAdjhEhEAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
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


   double beta_BMuPr;

   beta_BMuPr = Re(twoLoop*(-23.76*Power(g1,4)*MassB*MuPr - 15.36*Power(
      g1p,4)*MassBp*MuPr - 66*Power(g2,4)*MassWB*MuPr - 6*
      traceAdjgDTpYdconjYdTgD*MuPr - 6*traceAdjgDTpYuconjYuTgD*MuPr - 6*
      traceAdjYdTYdconjgDTpgD*MuPr - 6*traceAdjYuTYuconjgDTpgD*MuPr - 2*
      tracefuAdjhEThEAdjfu*MuPr - 36*tracegDAdjgDTgDAdjgD*MuPr - 6*
      tracegDAdjKappaTKappaAdjgD*MuPr - 2*tracehEAdjfuTfuAdjhE*MuPr - 12*
      tracehEAdjhEThEAdjhE*MuPr - 4*tracehEAdjhETYeAdjYe*MuPr - 2*
      tracehEAdjLambda12TLambda12AdjhE*MuPr - 6*traceKappaAdjgDTgDAdjKappa*MuPr
      - 2*traceLambda12AdjhEThEAdjLambda12*MuPr - 4*traceYeAdjYeThEAdjhE*MuPr
      - 18*traceAdjgDTgD*AbsSqr(SigmaL)*MuPr - 6*traceAdjhEThE*AbsSqr(SigmaL)*
      MuPr + 4*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax) + 2*AbsSqr(SigmaL))*BMuPhi*
      Conj(KappaPr)*SigmaL - 0.8*traceAdjgDTgD*MuPr*Sqr(g1) + 2.4*traceAdjhEThE
      *MuPr*Sqr(g1) + 0.8*MassB*tracegDAdjgD*MuPr*Sqr(g1) - 2.4*MassB*
      tracehEAdjhE*MuPr*Sqr(g1) - 9.6*MassB*AbsSqr(SigmaL)*MuPr*Sqr(g1) + 1.8*
      traceAdjgDTgD*MuPr*Sqr(g1p) + 0.6*traceAdjhEThE*MuPr*Sqr(g1p) - 1.8*
      MassBp*tracegDAdjgD*MuPr*Sqr(g1p) - 0.6*MassBp*tracehEAdjhE*MuPr*Sqr(g1p)
      - 6.4*MassBp*AbsSqr(SigmaL)*MuPr*Sqr(g1p) - 1.44*MassB*MuPr*Sqr(g1)*Sqr(
      g1p) - 1.44*MassBp*MuPr*Sqr(g1)*Sqr(g1p) - 48*MassWB*AbsSqr(SigmaL)*MuPr*
      Sqr(g2) - 3.6*MassB*MuPr*Sqr(g1)*Sqr(g2) - 3.6*MassWB*MuPr*Sqr(g1)*Sqr(g2
      ) - 2.4*MassBp*MuPr*Sqr(g1p)*Sqr(g2) - 2.4*MassWB*MuPr*Sqr(g1p)*Sqr(g2) +
      32*traceAdjgDTgD*MuPr*Sqr(g3) - 32*MassG*tracegDAdjgD*MuPr*Sqr(g3) -
      0.08*Power(g1p,4)*MassBp*MuPr*Sqr(QS) + 0.02*BMuPr*(297*Power(g1,4) + 192
      *Power(g1p,4) + 825*Power(g2,4) - 50*tracefuAdjhEhEAdjfu - 450*
      tracegDAdjgDgDAdjgD - 150*tracegDAdjgDTpYdconjYd - 150*
      tracegDAdjgDTpYuconjYu - 150*tracegDAdjKappaKappaAdjgD - 150*
      tracehEAdjhEhEAdjhE - 100*tracehEAdjhEYeAdjYe - 50*
      tracehEAdjLambda12Lambda12AdjhE + 60*tracehEAdjhE*Sqr(g1) + 15*
      tracehEAdjhE*Sqr(g1p) + 36*Sqr(g1)*Sqr(g1p) + 90*Sqr(g1)*Sqr(g2) + 60*Sqr
      (g1p)*Sqr(g2) + 10*AbsSqr(SigmaL)*(-75*tracegDAdjgD - 25*tracehEAdjhE -
      20*AbsSqr(KappaPr) - 10*AbsSqr(Sigmax) + 48*Sqr(g1) + 32*Sqr(g1p) + 240*
      Sqr(g2)) + tracegDAdjgD*(-20*Sqr(g1) + 45*Sqr(g1p) + 800*Sqr(g3)) + Power
      (g1p,4)*Sqr(QS) - 700*Sqr(Conj(SigmaL))*Sqr(SigmaL)) - 8*AbsSqr(SigmaL)*
      Conj(KappaPr)*MuPr*TKappaPr + 8*MuPhi*SigmaL*Sqr(Conj(KappaPr))*TKappaPr
      - 4*AbsSqr(SigmaL)*Conj(Sigmax)*MuPr*TSigmax + 4*MuPhi*Conj(KappaPr)*Conj
      (Sigmax)*SigmaL*TSigmax + 8*MuPhi*AbsSqr(SigmaL)*Conj(KappaPr)*TSigmaL -
      6*tracegDAdjgD*Conj(SigmaL)*MuPr*TSigmaL - 2*tracehEAdjhE*Conj(SigmaL)*
      MuPr*TSigmaL - 8*AbsSqr(KappaPr)*Conj(SigmaL)*MuPr*TSigmaL - 4*AbsSqr(
      Sigmax)*Conj(SigmaL)*MuPr*TSigmaL - 32*MuPr*SigmaL*Sqr(Conj(SigmaL))*
      TSigmaL));


   return beta_BMuPr;
}

} // namespace flexiblesusy
