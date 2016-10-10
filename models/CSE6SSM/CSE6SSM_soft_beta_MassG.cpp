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
 * Calculates the one-loop beta function of MassG.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassG_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassG;

   beta_MassG = Re(0);


   return beta_MassG;
}

/**
 * Calculates the two-loop beta function of MassG.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_MassG_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_MassG;

   beta_MassG = Re(2*twoLoop*Sqr(g3)*(4*traceAdjgDTgD + 2*
      traceAdjKappaTKappa + 4*traceAdjYdTYd + 4*traceAdjYuTYu - 4*MassG*
      tracegDAdjgD - 2*MassG*traceKappaAdjKappa - 4*MassG*traceYdAdjYd - 4*
      MassG*traceYuAdjYu + 3*MassB*Sqr(g1) + 3*MassG*Sqr(g1) + 3*MassBp*Sqr(g1p
      ) + 3*MassG*Sqr(g1p) + 9*MassG*Sqr(g2) + 9*MassWB*Sqr(g2) + 96*MassG*Sqr(
      g3)));


   return beta_MassG;
}

} // namespace flexiblesusy
