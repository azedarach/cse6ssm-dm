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

// File generated at Wed 3 Jun 2015 23:42:45

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gD.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_gD_one_loop(const Susy_traces& susy_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,3> beta_gD;

   beta_gD = (oneOver16PiSqr*(gD*(3*tracegDAdjgD + tracehEAdjhE + AbsSqr(
      SigmaL) - 0.4666666666666667*Sqr(g1) - 0.7*Sqr(g1p) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3)) + 3*(gD*gD.adjoint()*gD) + gD*(Kappa).adjoint(
      )*Kappa + Yd.transpose()*Yd.conjugate()*gD + Yu.transpose()*Yu.conjugate(
      )*gD)).real();


   return beta_gD;
}

/**
 * Calculates the two-loop beta function of gD.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_gD_two_loop(const Susy_traces& susy_traces) const
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
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_gD;

   beta_gD = (twoLoop*(gD*(4.588888888888889*Power(g1,4) + 6.825*Power(
      g1p,4) + 16.5*Power(g2,4) + 14.222222222222221*Power(g3,4) -
      tracefuAdjhEhEAdjfu - 9*tracegDAdjgDgDAdjgD - 3*tracegDAdjgDTpYdconjYd -
      3*tracegDAdjgDTpYuconjYu - 3*tracegDAdjKappaKappaAdjgD - 3*
      tracehEAdjhEhEAdjhE - 2*tracehEAdjhEYeAdjYe -
      tracehEAdjLambda12Lambda12AdjhE - 2*AbsSqr(KappaPr)*AbsSqr(SigmaL) -
      AbsSqr(Sigmax)*AbsSqr(SigmaL) - 0.4*tracegDAdjgD*Sqr(g1) + 1.2*
      tracehEAdjhE*Sqr(g1) + 0.9*tracegDAdjgD*Sqr(g1p) + 0.3*tracehEAdjhE*Sqr(
      g1p) + 0.6833333333333333*Sqr(g1)*Sqr(g1p) + Sqr(g1)*Sqr(g2) + 0.75*Sqr(
      g1p)*Sqr(g2) + 16*tracegDAdjgD*Sqr(g3) + 0.8888888888888888*Sqr(g1)*Sqr(
      g3) + 2.6666666666666665*Sqr(g1p)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) + 0.035*
      Power(g1p,4)*Sqr(QS) - 3*Sqr(Conj(SigmaL))*Sqr(SigmaL)) + (-9*
      tracegDAdjgD - 3*tracehEAdjhE - 3*AbsSqr(SigmaL) + 0.8*Sqr(g1) + 0.2*Sqr(
      g1p) + 6*Sqr(g2))*(gD*gD.adjoint()*gD) - 3*traceKappaAdjKappa*(gD*(Kappa)
      .adjoint()*Kappa) - 2*traceLambda12AdjLambda12*(gD*(Kappa).adjoint()*
      Kappa) - 2*AbsSqr(Lambdax)*(gD*(Kappa).adjoint()*Kappa) - AbsSqr(Sigmax)*
      (gD*(Kappa).adjoint()*Kappa) - 0.25*Sqr(g1p)*(gD*(Kappa).adjoint()*Kappa)
      + 0.05*Sqr(g1p)*Sqr(QS)*(gD*(Kappa).adjoint()*Kappa) - tracefdAdjfd*(
      Yd.transpose()*Yd.conjugate()*gD) - 3*traceYdAdjYd*(Yd.transpose()*
      Yd.conjugate()*gD) - traceYeAdjYe*(Yd.transpose()*Yd.conjugate()*gD) -
      AbsSqr(Lambdax)*(Yd.transpose()*Yd.conjugate()*gD) + 0.4*Sqr(g1)*(
      Yd.transpose()*Yd.conjugate()*gD) + 0.6*Sqr(g1p)*(Yd.transpose()*
      Yd.conjugate()*gD) - tracefuAdjfu*(Yu.transpose()*Yu.conjugate()*gD) - 3*
      traceYuAdjYu*(Yu.transpose()*Yu.conjugate()*gD) - AbsSqr(Lambdax)*(
      Yu.transpose()*Yu.conjugate()*gD) + 0.8*Sqr(g1)*(Yu.transpose()*
      Yu.conjugate()*gD) + 0.2*Sqr(g1p)*(Yu.transpose()*Yu.conjugate()*gD) - 4*
      (gD*gD.adjoint()*gD*gD.adjoint()*gD) - 2*(gD*gD.adjoint()*Yd.transpose()*
      Yd.conjugate()*gD) - 2*(gD*gD.adjoint()*Yu.transpose()*Yu.conjugate()*gD)
      - gD*(Kappa).adjoint()*Kappa*gD.adjoint()*gD - gD*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()*Kappa - 2*(Yd.transpose()*Yd.conjugate()*
      Yd.transpose()*Yd.conjugate()*gD) - 2*(Yu.transpose()*Yu.conjugate()*
      Yu.transpose()*Yu.conjugate()*gD))).real();


   return beta_gD;
}

} // namespace flexiblesusy
