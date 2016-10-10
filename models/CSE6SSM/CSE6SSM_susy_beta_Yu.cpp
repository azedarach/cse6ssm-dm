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
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(tracefuAdjfu + 3*traceYuAdjYu + AbsSqr(
      Lambdax) - 0.8666666666666667*Sqr(g1) - 0.3*Sqr(g1p) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3)) + Yu*Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu)
      + Yu*gD.conjugate()*gD.transpose())).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(Yu*(8.695555555555556*Power(g1,4) + 2.865*Power(
      g1p,4) + 16.5*Power(g2,4) + 14.222222222222221*Power(g3,4) - 2*
      tracefdAdjfdfuAdjfu - 3*tracefuAdjfufuAdjfu - tracefuAdjhEhEAdjfu -
      tracefuAdjLambda12Lambda12Adjfu - 3*tracegDAdjgDTpYuconjYu - 3*
      traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu + 0.8*traceYuAdjYu*Sqr(g1) +
      1.5*tracefuAdjfu*Sqr(g1p) - 0.3*traceYuAdjYu*Sqr(g1p) +
      0.5366666666666666*Sqr(g1)*Sqr(g1p) + Sqr(g1)*Sqr(g2) + 0.75*Sqr(g1p)*Sqr
      (g2) + 16*traceYuAdjYu*Sqr(g3) + 3.022222222222222*Sqr(g1)*Sqr(g3) +
      0.5333333333333333*Sqr(g1p)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) + 0.015*Power(g1p
      ,4)*Sqr(QS) + 0.05*AbsSqr(Lambdax)*(-20*tracefdAdjfd - 60*
      traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 60*traceYdAdjYd - 20*
      traceYeAdjYe - 20*AbsSqr(Sigmax) + 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS)) - 3*Sqr
      (Conj(Lambdax))*Sqr(Lambdax)) + 0.2*(-5*tracefdAdjfd - 15*traceYdAdjYd -
      5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(g1p))*(Yu*
      Yd.adjoint()*Yd) - 3*tracefuAdjfu*(Yu*Yu.adjoint()*Yu) - 9*traceYuAdjYu*(
      Yu*Yu.adjoint()*Yu) - 3*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*Yu) + 0.4*Sqr(g1
      )*(Yu*Yu.adjoint()*Yu) + 0.6*Sqr(g1p)*(Yu*Yu.adjoint()*Yu) + 6*Sqr(g2)*(
      Yu*Yu.adjoint()*Yu) - 3*tracegDAdjgD*(Yu*gD.conjugate()*gD.transpose()) -
      tracehEAdjhE*(Yu*gD.conjugate()*gD.transpose()) - AbsSqr(SigmaL)*(Yu*
      gD.conjugate()*gD.transpose()) + 0.4*Sqr(g1)*(Yu*gD.conjugate()*
      gD.transpose()) + 0.6*Sqr(g1p)*(Yu*gD.conjugate()*gD.transpose()) - 2*(Yu
      *Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu
      ) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 2*(Yu*gD.conjugate()*
      gD.transpose()*Yu.adjoint()*Yu) - 2*(Yu*gD.conjugate()*gD.transpose()*
      gD.conjugate()*gD.transpose()) - Yu*gD.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()*gD.transpose())).real();


   return beta_Yu;
}

} // namespace flexiblesusy
