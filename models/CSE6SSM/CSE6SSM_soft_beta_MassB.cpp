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
 * Calculates the one-loop beta function of MassB.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(19.2*MassB*oneOver16PiSqr*Sqr(g1));


   return beta_MassB;
}

/**
 * Calculates the two-loop beta function of MassB.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassB_two_loop(const Soft_traces& soft_traces) const
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


   double beta_MassB;

   beta_MassB = Re(0.08*twoLoop*Sqr(g1)*(30*traceAdjfdTfd + 30*
      traceAdjfuTfu + 70*traceAdjgDTgD + 90*traceAdjhEThE + 20*
      traceAdjKappaTKappa + 30*traceAdjLambda12TLambda12 + 70*traceAdjYdTYd +
      90*traceAdjYeTYe + 130*traceAdjYuTYu - 30*MassB*tracefdAdjfd - 30*MassB*
      tracefuAdjfu - 70*MassB*tracegDAdjgD - 90*MassB*tracehEAdjhE - 20*MassB*
      traceKappaAdjKappa - 30*MassB*traceLambda12AdjLambda12 - 70*MassB*
      traceYdAdjYd - 90*MassB*traceYeAdjYe - 130*MassB*traceYuAdjYu + 468*MassB
      *Sqr(g1) + 81*MassB*Sqr(g1p) + 81*MassBp*Sqr(g1p) + 270*MassB*Sqr(g2) +
      270*MassWB*Sqr(g2) + 600*MassB*Sqr(g3) + 600*MassG*Sqr(g3) - 30*Conj(
      Lambdax)*(MassB*Lambdax - TLambdax) - 30*Conj(SigmaL)*(MassB*SigmaL -
      TSigmaL)));


   return beta_MassB;
}

} // namespace flexiblesusy
