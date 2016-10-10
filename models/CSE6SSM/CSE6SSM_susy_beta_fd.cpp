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
 * Calculates the one-loop beta function of fd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_fd_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,2> beta_fd;

   beta_fd = (oneOver16PiSqr*(fd*(tracefdAdjfd + 3*traceYdAdjYd +
      traceYeAdjYe + AbsSqr(Lambdax) - 0.6*Sqr(g1) - 1.9*Sqr(g1p) - 3*Sqr(g2))
      + 3*(fd*fd.adjoint()*fd) + fd*Lambda12.conjugate()*(Lambda12).transpose()
      + 2*(fu*fu.adjoint()*fd))).real();


   return beta_fd;
}

/**
 * Calculates the two-loop beta function of fd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_susy_parameters::calc_beta_fd_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfdfdAdjfd = TRACE_STRUCT.tracefdAdjfdfdAdjfd;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;


   Eigen::Matrix<double,3,2> beta_fd;

   beta_fd = (twoLoop*(0.005*fd*(1188*Power(g1,4) + 3933*Power(g1p,4) +
      3300*Power(g2,4) - 600*tracefdAdjfdfdAdjfd - 400*tracefdAdjfdfuAdjfu -
      600*tracegDAdjgDTpYdconjYd - 400*tracehEAdjhEYeAdjYe - 200*
      traceLambda12AdjLambda12Tpfdconjfd - 1800*traceYdAdjYdYdAdjYd - 600*
      traceYdAdjYuYuAdjYd - 600*traceYeAdjYeYeAdjYe - 80*traceYdAdjYd*Sqr(g1) +
      240*traceYeAdjYe*Sqr(g1) + 200*tracefdAdjfd*Sqr(g1p) - 120*traceYdAdjYd*
      Sqr(g1p) - 40*traceYeAdjYe*Sqr(g1p) + 54*Sqr(g1)*Sqr(g1p) + 360*Sqr(g1)*
      Sqr(g2) + 390*Sqr(g1p)*Sqr(g2) + 3200*traceYdAdjYd*Sqr(g3) + 19*Power(g1p
      ,4)*Sqr(QS) + 10*AbsSqr(Lambdax)*(-20*tracefuAdjfu - 60*
      traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 60*traceYuAdjYu - 20*
      AbsSqr(Sigmax) - 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS)) - 600*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)) + (-3*tracefdAdjfd - 9*traceYdAdjYd - 3*traceYeAdjYe - 3*
      AbsSqr(Lambdax) + 1.2*Sqr(g1) + 0.3*Sqr(g1p) + 6*Sqr(g2))*(fd*fd.adjoint(
      )*fd) - 3*traceKappaAdjKappa*(fd*Lambda12.conjugate()*(Lambda12)
      .transpose()) - 2*traceLambda12AdjLambda12*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()) - 2*AbsSqr(Lambdax)*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()) - AbsSqr(Sigmax)*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()) + 0.25*Sqr(g1p)*(fd*Lambda12.conjugate()*(Lambda12
      ).transpose()) + 0.05*Sqr(g1p)*Sqr(QS)*(fd*Lambda12.conjugate()*(Lambda12
      ).transpose()) - 2*tracefuAdjfu*(fu*fu.adjoint()*fd) - 6*traceYuAdjYu*(fu
      *fu.adjoint()*fd) - 2*AbsSqr(Lambdax)*(fu*fu.adjoint()*fd) + 1.2*Sqr(g1)*
      (fu*fu.adjoint()*fd) - 1.2*Sqr(g1p)*(fu*fu.adjoint()*fd) + 6*Sqr(g2)*(fu*
      fu.adjoint()*fd) - 4*(fd*fd.adjoint()*fd*fd.adjoint()*fd) - 2*(fd*
      fd.adjoint()*fu*fu.adjoint()*fd) - fd*Lambda12.conjugate()*fu.transpose()
      *fu.conjugate()*(Lambda12).transpose() - fd*Lambda12.conjugate()*
      hE.transpose()*hE.conjugate()*(Lambda12).transpose() - 2*(fd*
      Lambda12.conjugate()*(Lambda12).transpose()*fd.adjoint()*fd) - fd*
      Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(
      Lambda12).transpose() - 2*(fu*fu.adjoint()*fu*fu.adjoint()*fd) - 2*(fu*
      hE.adjoint()*hE*fu.adjoint()*fd) - 2*(fu*(Lambda12).adjoint()*Lambda12*
      fu.adjoint()*fd))).real();


   return beta_fd;
}

} // namespace flexiblesusy
