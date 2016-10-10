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

// File generated at Wed 3 Jun 2015 23:42:43

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of hE.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_hE_one_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,2> beta_hE;

   beta_hE = (oneOver16PiSqr*(hE*(3*tracegDAdjgD + tracehEAdjhE + AbsSqr(
      SigmaL) - 1.8*Sqr(g1) - 0.7*Sqr(g1p) - 3*Sqr(g2)) + hE*fu.adjoint()*fu +
      3*(hE*hE.adjoint()*hE) + hE*(Lambda12).adjoint()*Lambda12 + 2*(Ye*
      Ye.adjoint()*hE))).real();


   return beta_hE;
}

/**
 * Calculates the two-loop beta function of hE.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_hE_two_loop(const Susy_traces& susy_traces) const
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
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,2> beta_hE;

   beta_hE = (twoLoop*(hE*(18.9*Power(g1,4) + 6.825*Power(g1p,4) + 16.5*
      Power(g2,4) - tracefuAdjhEhEAdjfu - 9*tracegDAdjgDgDAdjgD - 3*
      tracegDAdjgDTpYdconjYd - 3*tracegDAdjgDTpYuconjYu - 3*
      tracegDAdjKappaKappaAdjgD - 3*tracehEAdjhEhEAdjhE - 2*tracehEAdjhEYeAdjYe
      - tracehEAdjLambda12Lambda12AdjhE - 2*AbsSqr(KappaPr)*AbsSqr(SigmaL) -
      AbsSqr(Sigmax)*AbsSqr(SigmaL) - 0.4*tracegDAdjgD*Sqr(g1) + 1.2*
      tracehEAdjhE*Sqr(g1) + 0.9*tracegDAdjgD*Sqr(g1p) + 0.3*tracehEAdjhE*Sqr(
      g1p) + 0.15*Sqr(g1)*Sqr(g1p) + 1.8*Sqr(g1)*Sqr(g2) + 1.95*Sqr(g1p)*Sqr(g2
      ) + 16*tracegDAdjgD*Sqr(g3) + 0.035*Power(g1p,4)*Sqr(QS) - 3*Sqr(Conj(
      SigmaL))*Sqr(SigmaL)) + (-tracefuAdjfu - 3*traceYuAdjYu - AbsSqr(Lambdax)
      + Sqr(g1p))*(hE*fu.adjoint()*fu) - 9*tracegDAdjgD*(hE*hE.adjoint()*hE) -
      3*tracehEAdjhE*(hE*hE.adjoint()*hE) - 3*AbsSqr(SigmaL)*(hE*hE.adjoint()*
      hE) + Sqr(g1p)*(hE*hE.adjoint()*hE) + 6*Sqr(g2)*(hE*hE.adjoint()*hE) - 3*
      traceKappaAdjKappa*(hE*(Lambda12).adjoint()*Lambda12) - 2*
      traceLambda12AdjLambda12*(hE*(Lambda12).adjoint()*Lambda12) - 2*AbsSqr(
      Lambdax)*(hE*(Lambda12).adjoint()*Lambda12) - AbsSqr(Sigmax)*(hE*(
      Lambda12).adjoint()*Lambda12) - 0.25*Sqr(g1p)*(hE*(Lambda12).adjoint()*
      Lambda12) + 0.05*Sqr(g1p)*Sqr(QS)*(hE*(Lambda12).adjoint()*Lambda12) - 2*
      tracefdAdjfd*(Ye*Ye.adjoint()*hE) - 6*traceYdAdjYd*(Ye*Ye.adjoint()*hE) -
      2*traceYeAdjYe*(Ye*Ye.adjoint()*hE) - 2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()
      *hE) - 1.2*Sqr(g1)*(Ye*Ye.adjoint()*hE) + 1.2*Sqr(g1p)*(Ye*Ye.adjoint()*
      hE) + 6*Sqr(g2)*(Ye*Ye.adjoint()*hE) - 2*(hE*fu.adjoint()*fd*fd.adjoint()
      *fu) - 2*(hE*fu.adjoint()*fu*fu.adjoint()*fu) - 2*(hE*fu.adjoint()*fu*
      hE.adjoint()*hE) - 4*(hE*hE.adjoint()*hE*hE.adjoint()*hE) - 2*(hE*
      hE.adjoint()*Ye*Ye.adjoint()*hE) - 2*(hE*(Lambda12).adjoint()*Lambda12*
      hE.adjoint()*hE) - hE*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*
      Lambda12 - hE*(Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12
      - 2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*hE))).real();


   return beta_hE;
}

} // namespace flexiblesusy
