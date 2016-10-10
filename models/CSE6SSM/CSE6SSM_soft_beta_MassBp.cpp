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

// File generated at Wed 3 Jun 2015 23:43:44

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassBp.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassBp_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassBp;

   beta_MassBp = Re(0.1*MassBp*oneOver16PiSqr*Sqr(g1p)*(188 + Sqr(QS)));


   return beta_MassBp;
}

/**
 * Calculates the two-loop beta function of MassBp.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassBp_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_MassBp;

   beta_MassBp = Re(0.02*twoLoop*Sqr(g1p)*(380*traceAdjfdTfd + 380*
      traceAdjfuTfu + 420*traceAdjgDTgD + 140*traceAdjhEThE + 195*
      traceAdjKappaTKappa + 130*traceAdjLambda12TLambda12 + 420*traceAdjYdTYd +
      140*traceAdjYeTYe + 180*traceAdjYuTYu - 380*MassBp*tracefdAdjfd - 380*
      MassBp*tracefuAdjfu - 420*MassBp*tracegDAdjgD - 140*MassBp*tracehEAdjhE -
      195*MassBp*traceKappaAdjKappa - 130*MassBp*traceLambda12AdjLambda12 -
      420*MassBp*traceYdAdjYd - 140*MassBp*traceYeAdjYe - 180*MassBp*
      traceYuAdjYu - 80*MassBp*AbsSqr(SigmaL) + 324*MassB*Sqr(g1) + 324*MassBp*
      Sqr(g1) + 1832*MassBp*Sqr(g1p) + MassBp*Power(QS,4)*Sqr(g1p) + 1020*
      MassBp*Sqr(g2) + 1020*MassWB*Sqr(g2) + 2400*MassBp*Sqr(g3) + 2400*MassG*
      Sqr(g3) + 15*traceAdjKappaTKappa*Sqr(QS) + 10*traceAdjLambda12TLambda12*
      Sqr(QS) - 15*MassBp*traceKappaAdjKappa*Sqr(QS) - 10*MassBp*
      traceLambda12AdjLambda12*Sqr(QS) - 10*Conj(Lambdax)*(13 + Sqr(QS))*(
      MassBp*Lambdax - TLambdax) - 10*Conj(Sigmax)*Sqr(QS)*(MassBp*Sigmax -
      TSigmax) + 80*Conj(SigmaL)*TSigmaL));


   return beta_MassBp;
}

} // namespace flexiblesusy
