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

// File generated at Wed 3 Jun 2015 23:42:49

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vsb.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vsb_one_loop(const Susy_traces& susy_traces) const
{


   double beta_vsb;

   beta_vsb = Re(oneOver16PiSqr*(-(vsb*AbsSqr(Sigmax)) + 0.05*vsb*Sqr(g1p
      )*Sqr(QS)));


   return beta_vsb;
}

/**
 * Calculates the two-loop beta function of vsb.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_vsb_two_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_vsb;

   beta_vsb = Re(-0.0025*twoLoop*vsb*(Power(g1p,4)*Sqr(QS)*(94 + Sqr(QS))
      + 20*AbsSqr(Sigmax)*(-60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 40*AbsSqr(KappaPr) - 40*AbsSqr(Lambdax) - 40*
      AbsSqr(SigmaL) + Sqr(g1p)*Sqr(QS)) - 800*Sqr(Conj(Sigmax))*Sqr(Sigmax)));


   return beta_vsb;
}

} // namespace flexiblesusy
