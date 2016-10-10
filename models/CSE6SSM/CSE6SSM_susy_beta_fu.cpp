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

// File generated at Wed 3 Jun 2015 23:42:47

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of fu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_fu_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,2> beta_fu;

   beta_fu = (oneOver16PiSqr*(fu*(tracefuAdjfu + 3*traceYuAdjYu + AbsSqr(
      Lambdax) - 0.6*Sqr(g1) - 1.9*Sqr(g1p) - 3*Sqr(g2)) + 2*(fd*fd.adjoint()*
      fu) + 3*(fu*fu.adjoint()*fu) + fu*hE.adjoint()*hE + fu*(Lambda12).adjoint
      ()*Lambda12)).real();


   return beta_fu;
}

/**
 * Calculates the two-loop beta function of fu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_fu_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,2> beta_fu;

   beta_fu = (twoLoop*(0.005*fu*(1188*Power(g1,4) + 3933*Power(g1p,4) +
      3300*Power(g2,4) - 400*tracefdAdjfdfuAdjfu - 600*tracefuAdjfufuAdjfu -
      200*tracefuAdjhEhEAdjfu - 200*tracefuAdjLambda12Lambda12Adjfu - 600*
      tracegDAdjgDTpYuconjYu - 600*traceYdAdjYuYuAdjYd - 1800*
      traceYuAdjYuYuAdjYu + 160*traceYuAdjYu*Sqr(g1) + 300*tracefuAdjfu*Sqr(g1p
      ) - 60*traceYuAdjYu*Sqr(g1p) + 54*Sqr(g1)*Sqr(g1p) + 360*Sqr(g1)*Sqr(g2)
      + 390*Sqr(g1p)*Sqr(g2) + 3200*traceYuAdjYu*Sqr(g3) + 19*Power(g1p,4)*Sqr(
      QS) + 10*AbsSqr(Lambdax)*(-20*tracefdAdjfd - 60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 60*traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(
      Sigmax) + 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS)) - 600*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)) + (-2*tracefdAdjfd - 6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr
      (Lambdax) + 1.2*Sqr(g1) - 1.2*Sqr(g1p) + 6*Sqr(g2))*(fd*fd.adjoint()*fu)
      - 3*tracefuAdjfu*(fu*fu.adjoint()*fu) - 9*traceYuAdjYu*(fu*fu.adjoint()*
      fu) - 3*AbsSqr(Lambdax)*(fu*fu.adjoint()*fu) + 1.2*Sqr(g1)*(fu*fu.adjoint
      ()*fu) - 0.2*Sqr(g1p)*(fu*fu.adjoint()*fu) + 6*Sqr(g2)*(fu*fu.adjoint()*
      fu) - 3*tracegDAdjgD*(fu*hE.adjoint()*hE) - tracehEAdjhE*(fu*hE.adjoint()
      *hE) - AbsSqr(SigmaL)*(fu*hE.adjoint()*hE) + 1.2*Sqr(g1)*(fu*hE.adjoint()
      *hE) - 0.2*Sqr(g1p)*(fu*hE.adjoint()*hE) - 3*traceKappaAdjKappa*(fu*(
      Lambda12).adjoint()*Lambda12) - 2*traceLambda12AdjLambda12*(fu*(Lambda12)
      .adjoint()*Lambda12) - 2*AbsSqr(Lambdax)*(fu*(Lambda12).adjoint()*
      Lambda12) - AbsSqr(Sigmax)*(fu*(Lambda12).adjoint()*Lambda12) - 0.25*Sqr(
      g1p)*(fu*(Lambda12).adjoint()*Lambda12) + 0.05*Sqr(g1p)*Sqr(QS)*(fu*(
      Lambda12).adjoint()*Lambda12) - 2*(fd*fd.adjoint()*fd*fd.adjoint()*fu) -
      2*(fd*Lambda12.conjugate()*(Lambda12).transpose()*fd.adjoint()*fu) - 2*(
      fu*fu.adjoint()*fd*fd.adjoint()*fu) - 4*(fu*fu.adjoint()*fu*fu.adjoint()*
      fu) - 2*(fu*hE.adjoint()*hE*fu.adjoint()*fu) - 2*(fu*hE.adjoint()*hE*
      hE.adjoint()*hE) - 2*(fu*hE.adjoint()*Ye*Ye.adjoint()*hE) - 2*(fu*(
      Lambda12).adjoint()*Lambda12*fu.adjoint()*fu) - fu*(Lambda12).adjoint()*
      Lambda12*(Lambda12).adjoint()*Lambda12 - fu*(Lambda12).adjoint()*
      fd.transpose()*fd.conjugate()*Lambda12)).real();


   return beta_fu;
}

} // namespace flexiblesusy
