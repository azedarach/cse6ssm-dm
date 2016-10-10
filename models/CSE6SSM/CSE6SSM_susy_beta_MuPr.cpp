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
 * Calculates the one-loop beta function of MuPr.
 *
 * @return one-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_MuPr_one_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   double beta_MuPr;

   beta_MuPr = Re(0.2*oneOver16PiSqr*MuPr*(15*tracegDAdjgD + 5*
      tracehEAdjhE + 10*AbsSqr(SigmaL) - 3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2)));


   return beta_MuPr;
}

/**
 * Calculates the two-loop beta function of MuPr.
 *
 * @return two-loop beta function
 */
double CSE6SSM_susy_parameters::calc_beta_MuPr_two_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
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


   double beta_MuPr;

   beta_MuPr = Re(0.02*twoLoop*MuPr*(297*Power(g1,4) + 192*Power(g1p,4) +
      825*Power(g2,4) - 50*tracefuAdjhEhEAdjfu - 450*tracegDAdjgDgDAdjgD - 150
      *tracegDAdjgDTpYdconjYd - 150*tracegDAdjgDTpYuconjYu - 150*
      tracegDAdjKappaKappaAdjgD - 150*tracehEAdjhEhEAdjhE - 100*
      tracehEAdjhEYeAdjYe - 50*tracehEAdjLambda12Lambda12AdjhE - 150*
      tracegDAdjgD*AbsSqr(SigmaL) - 50*tracehEAdjhE*AbsSqr(SigmaL) - 200*AbsSqr
      (KappaPr)*AbsSqr(SigmaL) - 100*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 20*
      tracegDAdjgD*Sqr(g1) + 60*tracehEAdjhE*Sqr(g1) + 45*tracegDAdjgD*Sqr(g1p)
      + 15*tracehEAdjhE*Sqr(g1p) + 36*Sqr(g1)*Sqr(g1p) + 90*Sqr(g1)*Sqr(g2) +
      60*Sqr(g1p)*Sqr(g2) + 800*tracegDAdjgD*Sqr(g3) + Power(g1p,4)*Sqr(QS) -
      300*Sqr(Conj(SigmaL))*Sqr(SigmaL)));


   return beta_MuPr;
}

} // namespace flexiblesusy
