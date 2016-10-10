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
 * Calculates the one-loop beta function of BMuPhi.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_BMuPhi_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPhi;

   beta_BMuPhi = Re(oneOver16PiSqr*(2*(4*AbsSqr(KappaPr) + AbsSqr(Sigmax)
      + 2*AbsSqr(SigmaL))*BMuPhi - 8*BMuPr*Conj(SigmaL)*KappaPr + 4*MuPhi*(2*
      Conj(KappaPr)*TKappaPr + Conj(Sigmax)*TSigmax + 2*Conj(SigmaL)*TSigmaL)))
      ;


   return beta_BMuPhi;
}

/**
 * Calculates the two-loop beta function of BMuPhi.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_BMuPhi_two_loop(const Soft_traces& soft_traces) const
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


   double beta_BMuPhi;

   beta_BMuPhi = Re(twoLoop*(-0.2*BMuPhi*(80*AbsSqr(KappaPr)*(AbsSqr(
      Sigmax) + 2*AbsSqr(SigmaL)) + 4*AbsSqr(SigmaL)*(15*tracegDAdjgD + 5*
      tracehEAdjhE + 10*AbsSqr(SigmaL) - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2)) +
      AbsSqr(Sigmax)*(30*traceKappaAdjKappa + 20*traceLambda12AdjLambda12 + 20
      *AbsSqr(Lambdax) - Sqr(g1p)*Sqr(QS)) + 160*Sqr(Conj(KappaPr))*Sqr(KappaPr
      ) + 20*Sqr(Conj(Sigmax))*Sqr(Sigmax)) - 0.4*(40*MuPhi*Sigmax*Sqr(Conj(
      Sigmax))*TSigmax + MuPhi*Conj(Sigmax)*(30*traceAdjKappaTKappa*Sigmax + 20
      *traceAdjLambda12TLambda12*Sigmax + MassBp*Sigmax*Sqr(g1p)*Sqr(QS) + 30*
      traceKappaAdjKappa*TSigmax + 20*traceLambda12AdjLambda12*TSigmax - Sqr(
      g1p)*Sqr(QS)*TSigmax + 20*Conj(KappaPr)*(2*Sigmax*TKappaPr + 3*KappaPr*
      TSigmax) + 20*Conj(Lambdax)*(Sigmax*TLambdax + Lambdax*TSigmax)) - 4*(-50
      *MuPhi*KappaPr*Sqr(Conj(KappaPr))*TKappaPr + 10*Sqr(Conj(SigmaL))*(BMuPr*
      KappaPr*SigmaL + (KappaPr*MuPr - 2*MuPhi*SigmaL)*TSigmaL) + Conj(SigmaL)*
      (15*traceAdjgDTgD*KappaPr*MuPr + 5*traceAdjhEThE*KappaPr*MuPr - 15*MuPhi*
      traceAdjgDTgD*SigmaL - 5*MuPhi*traceAdjhEThE*SigmaL + 12*MassB*KappaPr*
      MuPr*Sqr(g1) - 3*MassB*MuPhi*SigmaL*Sqr(g1) + 8*MassBp*KappaPr*MuPr*Sqr(
      g1p) - 2*MassBp*MuPhi*SigmaL*Sqr(g1p) + BMuPr*KappaPr*(15*tracegDAdjgD +
      5*tracehEAdjhE - 12*Sqr(g1) - 8*Sqr(g1p) - 60*Sqr(g2)) + 60*MassWB*
      KappaPr*MuPr*Sqr(g2) - 15*MassWB*MuPhi*SigmaL*Sqr(g2) - 15*MuPhi*
      tracegDAdjgD*TSigmaL - 5*MuPhi*tracehEAdjhE*TSigmaL + 3*MuPhi*Sqr(g1)*
      TSigmaL + 2*MuPhi*Sqr(g1p)*TSigmaL + 15*MuPhi*Sqr(g2)*TSigmaL - 10*MuPhi*
      Conj(KappaPr)*(2*SigmaL*TKappaPr + 3*KappaPr*TSigmaL))))));


   return beta_BMuPhi;
}

} // namespace flexiblesusy
