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
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(tracefdAdjfd + 3*traceYdAdjYd +
      traceYeAdjYe + AbsSqr(Lambdax) - 1.8*Sqr(g1) - 0.7*Sqr(g1p) - 3*Sqr(g2))
      + 2*(hE*hE.adjoint()*Ye) + 3*(Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
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
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(Ye*(18.9*Power(g1,4) + 6.825*Power(g1p,4) + 16.5*
      Power(g2,4) - 3*tracefdAdjfdfdAdjfd - 2*tracefdAdjfdfuAdjfu - 3*
      tracegDAdjgDTpYdconjYd - 2*tracehEAdjhEYeAdjYe -
      traceLambda12AdjLambda12Tpfdconjfd - 9*traceYdAdjYdYdAdjYd - 3*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - 0.4*traceYdAdjYd*Sqr(g1) +
      1.2*traceYeAdjYe*Sqr(g1) + tracefdAdjfd*Sqr(g1p) - 0.6*traceYdAdjYd*Sqr(
      g1p) - 0.2*traceYeAdjYe*Sqr(g1p) + 0.15*Sqr(g1)*Sqr(g1p) + 1.8*Sqr(g1)*
      Sqr(g2) + 1.95*Sqr(g1p)*Sqr(g2) + 16*traceYdAdjYd*Sqr(g3) + 0.035*Power(
      g1p,4)*Sqr(QS) + 0.05*AbsSqr(Lambdax)*(-20*tracefuAdjfu - 60*
      traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 60*traceYuAdjYu - 20*
      AbsSqr(Sigmax) - 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS)) - 3*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)) - 0.4*(15*tracegDAdjgD + 5*tracehEAdjhE + 5*AbsSqr(SigmaL)
      + 3*Sqr(g1) - 3*Sqr(g1p) - 15*Sqr(g2))*(hE*hE.adjoint()*Ye) - 3*
      tracefdAdjfd*(Ye*Ye.adjoint()*Ye) - 9*traceYdAdjYd*(Ye*Ye.adjoint()*Ye) -
      3*traceYeAdjYe*(Ye*Ye.adjoint()*Ye) - 3*AbsSqr(Lambdax)*(Ye*Ye.adjoint()
      *Ye) + 1.5*Sqr(g1p)*(Ye*Ye.adjoint()*Ye) + 6*Sqr(g2)*(Ye*Ye.adjoint()*Ye)
      - 2*(hE*fu.adjoint()*fu*hE.adjoint()*Ye) - 2*(hE*hE.adjoint()*hE*
      hE.adjoint()*Ye) - 2*(hE*(Lambda12).adjoint()*Lambda12*hE.adjoint()*Ye) -
      2*(Ye*Ye.adjoint()*hE*hE.adjoint()*Ye) - 4*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

} // namespace flexiblesusy
