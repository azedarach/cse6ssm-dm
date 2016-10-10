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

// File generated at Wed 3 Jun 2015 23:42:48

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MuPhi.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_MuPhi_one_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPhi;

   beta_MuPhi = Re(2*MuPhi*oneOver16PiSqr*(2*AbsSqr(KappaPr) + AbsSqr(
      Sigmax) + 2*AbsSqr(SigmaL)));


   return beta_MuPhi;
}

/**
 * Calculates the two-loop beta function of MuPhi.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_MuPhi_two_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_MuPhi;

   beta_MuPhi = Re(-0.2*MuPhi*twoLoop*(40*AbsSqr(KappaPr)*(AbsSqr(Sigmax)
      + 2*AbsSqr(SigmaL)) + 4*AbsSqr(SigmaL)*(15*tracegDAdjgD + 5*tracehEAdjhE
      + 10*AbsSqr(SigmaL) - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2)) + AbsSqr(
      Sigmax)*(30*traceKappaAdjKappa + 20*traceLambda12AdjLambda12 + 20*AbsSqr(
      Lambdax) - Sqr(g1p)*Sqr(QS)) + 80*Sqr(Conj(KappaPr))*Sqr(KappaPr) + 20*
      Sqr(Conj(Sigmax))*Sqr(Sigmax)));


   return beta_MuPhi;
}

} // namespace flexiblesusy
