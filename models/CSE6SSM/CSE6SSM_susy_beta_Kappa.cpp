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
 * Calculates the one-loop beta function of Kappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Kappa_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (oneOver16PiSqr*(Kappa*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) + AbsSqr(Sigmax) -
      0.26666666666666666*Sqr(g1) - 0.65*Sqr(g1p) - 5.333333333333333*Sqr(g3) -
      0.05*Sqr(g1p)*Sqr(QS)) + 2*(Kappa*gD.adjoint()*gD + Kappa*(Kappa)
      .adjoint()*Kappa))).real();


   return beta_Kappa;
}

/**
 * Calculates the two-loop beta function of Kappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_susy_parameters::calc_beta_Kappa_two_loop(const Susy_traces& susy_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12Tpfdconjfd =
      TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (twoLoop*(Kappa*(2.5955555555555554*Power(g1,4) + 6.3525*
      Power(g1p,4) + 14.222222222222221*Power(g3,4) + 0.005*Power(g1p,4)*Power(
      QS,4) - 2*tracefuAdjLambda12Lambda12Adjfu - 6*tracegDAdjKappaKappaAdjgD -
      2*tracehEAdjLambda12Lambda12AdjhE - 6*traceKappaAdjKappaKappaAdjKappa -
      4*traceLambda12AdjLambda12Lambda12AdjLambda12 - 2*
      traceLambda12AdjLambda12Tpfdconjfd - 2*AbsSqr(KappaPr)*AbsSqr(Sigmax) - 2
      *AbsSqr(Sigmax)*AbsSqr(SigmaL) + 0.8*traceKappaAdjKappa*Sqr(g1) + 1.2*
      traceLambda12AdjLambda12*Sqr(g1) + 1.95*traceKappaAdjKappa*Sqr(g1p) + 1.3
      *traceLambda12AdjLambda12*Sqr(g1p) + 0.25333333333333335*Sqr(g1)*Sqr(g1p)
      + 6*traceLambda12AdjLambda12*Sqr(g2) + 16*traceKappaAdjKappa*Sqr(g3) +
      1.4222222222222223*Sqr(g1)*Sqr(g3) + 3.466666666666667*Sqr(g1p)*Sqr(g3) +
      0.5025*Power(g1p,4)*Sqr(QS) - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) -
      0.1*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS) + 0.1*AbsSqr(Lambdax)*(-20
      *tracefdAdjfd - 20*tracefuAdjfu - 60*traceYdAdjYd - 20*traceYeAdjYe - 60*
      traceYuAdjYu + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2) - Sqr(g1p)*Sqr(QS))
      - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)) +
      0.4*(-15*tracegDAdjgD - 5*tracehEAdjhE - 5*AbsSqr(SigmaL) + Sqr(g1) - Sqr
      (g1p) + 15*Sqr(g2))*(Kappa*gD.adjoint()*gD) - 6*traceKappaAdjKappa*(Kappa
      *(Kappa).adjoint()*Kappa) - 4*traceLambda12AdjLambda12*(Kappa*(Kappa)
      .adjoint()*Kappa) - 4*AbsSqr(Lambdax)*(Kappa*(Kappa).adjoint()*Kappa) - 2
      *AbsSqr(Sigmax)*(Kappa*(Kappa).adjoint()*Kappa) + 0.1*Sqr(g1p)*Sqr(QS)*(
      Kappa*(Kappa).adjoint()*Kappa) - 2*(Kappa*gD.adjoint()*gD*gD.adjoint()*gD
      ) - 2*(Kappa*gD.adjoint()*gD*(Kappa).adjoint()*Kappa) - 2*(Kappa*
      gD.adjoint()*Yd.transpose()*Yd.conjugate()*gD) - 2*(Kappa*gD.adjoint()*
      Yu.transpose()*Yu.conjugate()*gD) - 2*(Kappa*(Kappa).adjoint()*Kappa*(
      Kappa).adjoint()*Kappa))).real();


   return beta_Kappa;
}

} // namespace flexiblesusy
