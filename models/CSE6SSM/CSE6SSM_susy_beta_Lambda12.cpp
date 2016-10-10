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

// File generated at Wed 3 Jun 2015 23:42:46

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_susy_parameters::calc_beta_Lambda12_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = (oneOver16PiSqr*(Lambda12*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) + AbsSqr(Sigmax) - 0.6*Sqr(
      g1) - 0.65*Sqr(g1p) - 3*Sqr(g2) - 0.05*Sqr(g1p)*Sqr(QS)) + Lambda12*
      fu.adjoint()*fu + Lambda12*hE.adjoint()*hE + 2*(Lambda12*(Lambda12)
      .adjoint()*Lambda12) + fd.transpose()*fd.conjugate()*Lambda12)).real();


   return beta_Lambda12;
}

/**
 * Calculates the two-loop beta function of Lambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_susy_parameters::calc_beta_Lambda12_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = (twoLoop*(0.0025*Lambda12*(2376*Power(g1,4) + 2541*
      Power(g1p,4) + 6600*Power(g2,4) + 2*Power(g1p,4)*Power(QS,4) - 800*
      tracefuAdjLambda12Lambda12Adjfu - 2400*tracegDAdjKappaKappaAdjgD - 800*
      tracehEAdjLambda12Lambda12AdjhE - 2400*traceKappaAdjKappaKappaAdjKappa -
      1600*traceLambda12AdjLambda12Lambda12AdjLambda12 - 800*
      traceLambda12AdjLambda12Tpfdconjfd - 800*AbsSqr(KappaPr)*AbsSqr(Sigmax) -
      800*AbsSqr(Sigmax)*AbsSqr(SigmaL) + 320*traceKappaAdjKappa*Sqr(g1) + 480
      *traceLambda12AdjLambda12*Sqr(g1) + 780*traceKappaAdjKappa*Sqr(g1p) + 520
      *traceLambda12AdjLambda12*Sqr(g1p) + 108*Sqr(g1)*Sqr(g1p) + 2400*
      traceLambda12AdjLambda12*Sqr(g2) + 720*Sqr(g1)*Sqr(g2) + 780*Sqr(g1p)*Sqr
      (g2) + 6400*traceKappaAdjKappa*Sqr(g3) + 201*Power(g1p,4)*Sqr(QS) - 60*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) - 40*traceLambda12AdjLambda12*Sqr(g1p
      )*Sqr(QS) + 40*AbsSqr(Lambdax)*(-20*tracefdAdjfd - 20*tracefuAdjfu - 60*
      traceYdAdjYd - 20*traceYeAdjYe - 60*traceYuAdjYu + 12*Sqr(g1) + 13*Sqr(
      g1p) + 60*Sqr(g2) - Sqr(g1p)*Sqr(QS)) - 1600*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 800*Sqr(Conj(Sigmax))*Sqr(Sigmax)) + (-tracefuAdjfu - 3*
      traceYuAdjYu - AbsSqr(Lambdax) + Sqr(g1p))*(Lambda12*fu.adjoint()*fu) - 3
      *tracegDAdjgD*(Lambda12*hE.adjoint()*hE) - tracehEAdjhE*(Lambda12*
      hE.adjoint()*hE) - AbsSqr(SigmaL)*(Lambda12*hE.adjoint()*hE) + 1.2*Sqr(g1
      )*(Lambda12*hE.adjoint()*hE) - 0.2*Sqr(g1p)*(Lambda12*hE.adjoint()*hE) -
      6*traceKappaAdjKappa*(Lambda12*(Lambda12).adjoint()*Lambda12) - 4*
      traceLambda12AdjLambda12*(Lambda12*(Lambda12).adjoint()*Lambda12) - 4*
      AbsSqr(Lambdax)*(Lambda12*(Lambda12).adjoint()*Lambda12) - 2*AbsSqr(
      Sigmax)*(Lambda12*(Lambda12).adjoint()*Lambda12) + 0.1*Sqr(g1p)*Sqr(QS)*(
      Lambda12*(Lambda12).adjoint()*Lambda12) - tracefdAdjfd*(fd.transpose()*
      fd.conjugate()*Lambda12) - 3*traceYdAdjYd*(fd.transpose()*fd.conjugate()*
      Lambda12) - traceYeAdjYe*(fd.transpose()*fd.conjugate()*Lambda12) -
      AbsSqr(Lambdax)*(fd.transpose()*fd.conjugate()*Lambda12) + 1.5*Sqr(g1p)*(
      fd.transpose()*fd.conjugate()*Lambda12) - 2*(Lambda12*fu.adjoint()*fd*
      fd.adjoint()*fu) - 2*(Lambda12*fu.adjoint()*fu*fu.adjoint()*fu) -
      Lambda12*fu.adjoint()*fu*(Lambda12).adjoint()*Lambda12 - 2*(Lambda12*
      hE.adjoint()*hE*hE.adjoint()*hE) - Lambda12*hE.adjoint()*hE*(Lambda12)
      .adjoint()*Lambda12 - 2*(Lambda12*hE.adjoint()*Ye*Ye.adjoint()*hE) - 2*(
      Lambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12) -
      Lambda12*(Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12 - 2*
      (fd.transpose()*fd.conjugate()*fd.transpose()*fd.conjugate()*Lambda12) -
      2*(fd.transpose()*fu.conjugate()*fu.transpose()*fd.conjugate()*Lambda12))
      ).real();


   return beta_Lambda12;
}

} // namespace flexiblesusy
