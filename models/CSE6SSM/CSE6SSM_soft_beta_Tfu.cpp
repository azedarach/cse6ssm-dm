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

// File generated at Wed 3 Jun 2015 23:43:11

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of Tfu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_Tfu_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,2> beta_Tfu;

   beta_Tfu = (oneOver16PiSqr*(tracefuAdjfu*Tfu + 3*traceYuAdjYu*Tfu +
      AbsSqr(Lambdax)*Tfu - 0.6*Sqr(g1)*Tfu - 1.9*Sqr(g1p)*Tfu - 3*Sqr(g2)*Tfu
      + fu*(2*traceAdjfuTfu + 6*traceAdjYuTYu + 1.2*MassB*Sqr(g1) + 3.8*MassBp*
      Sqr(g1p) + 6*MassWB*Sqr(g2) + 2*Conj(Lambdax)*TLambdax) + 2*(fd*
      fd.adjoint()*Tfu) + 4*(fu*fu.adjoint()*Tfu) + 2*(fu*hE.adjoint()*ThE) + 2
      *(fu*(Lambda12).adjoint()*TLambda12) + 4*(Tfd*fd.adjoint()*fu) + 5*(Tfu*
      fu.adjoint()*fu) + Tfu*hE.adjoint()*hE + Tfu*(Lambda12).adjoint()*
      Lambda12)).real();


   return beta_Tfu;
}

/**
 * Calculates the two-loop beta function of Tfu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_Tfu_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double tracefdAdjfdTfuAdjfu = TRACE_STRUCT.tracefdAdjfdTfuAdjfu;
   const double tracefuAdjfuTfdAdjfd = TRACE_STRUCT.tracefuAdjfuTfdAdjfd;
   const double tracefuAdjfuTfuAdjfu = TRACE_STRUCT.tracefuAdjfuTfuAdjfu;
   const double tracefuAdjhEThEAdjfu = TRACE_STRUCT.tracefuAdjhEThEAdjfu;
   const double tracefuAdjLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12Adjfu;
   const double tracehEAdjfuTfuAdjhE = TRACE_STRUCT.tracehEAdjfuTfuAdjhE;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceLambda12AdjfuTfuAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjfuTfuAdjLambda12;
   const double traceAdjgDTpYuconjYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjYuTgD;
   const double traceAdjYuTYuconjgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjgDTpgD;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double tracefdAdjfdfuAdjfu = TRACE_STRUCT.tracefdAdjfdfuAdjfu;
   const double tracefuAdjfufuAdjfu = TRACE_STRUCT.tracefuAdjfufuAdjfu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracefuAdjLambda12Lambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,2> beta_Tfu;

   beta_Tfu = (twoLoop*(5.94*Power(g1,4)*Tfu + 19.665*Power(g1p,4)*Tfu +
      16.5*Power(g2,4)*Tfu - 2*tracefdAdjfdfuAdjfu*Tfu - 3*tracefuAdjfufuAdjfu*
      Tfu - tracefuAdjhEhEAdjfu*Tfu - tracefuAdjLambda12Lambda12Adjfu*Tfu - 3*
      tracegDAdjgDTpYuconjYu*Tfu - 3*traceYdAdjYuYuAdjYd*Tfu - 9*
      traceYuAdjYuYuAdjYu*Tfu - tracefdAdjfd*AbsSqr(Lambdax)*Tfu - 3*
      traceKappaAdjKappa*AbsSqr(Lambdax)*Tfu - 2*traceLambda12AdjLambda12*
      AbsSqr(Lambdax)*Tfu - 3*traceYdAdjYd*AbsSqr(Lambdax)*Tfu - traceYeAdjYe*
      AbsSqr(Lambdax)*Tfu - AbsSqr(Lambdax)*AbsSqr(Sigmax)*Tfu + 0.8*
      traceYuAdjYu*Sqr(g1)*Tfu + 1.5*tracefuAdjfu*Sqr(g1p)*Tfu - 0.3*
      traceYuAdjYu*Sqr(g1p)*Tfu + 0.25*AbsSqr(Lambdax)*Sqr(g1p)*Tfu + 0.27*Sqr(
      g1)*Sqr(g1p)*Tfu + 1.8*Sqr(g1)*Sqr(g2)*Tfu + 1.95*Sqr(g1p)*Sqr(g2)*Tfu +
      16*traceYuAdjYu*Sqr(g3)*Tfu + 0.095*Power(g1p,4)*Sqr(QS)*Tfu + 0.05*
      AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*Tfu - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*
      Tfu - 0.02*fu*(1188*Power(g1,4)*MassB + 3933*Power(g1p,4)*MassBp + 3300*
      Power(g2,4)*MassWB + 300*traceAdjgDTpYuconjYuTgD + 300*
      traceAdjYuTYuconjgDTpgD + 200*tracefdAdjfdTfuAdjfu + 200*
      tracefuAdjfuTfdAdjfd + 600*tracefuAdjfuTfuAdjfu + 100*
      tracefuAdjhEThEAdjfu + 100*tracefuAdjLambda12TLambda12Adjfu + 100*
      tracehEAdjfuTfuAdjhE + 100*traceLambda12AdjfuTfuAdjLambda12 + 300*
      traceYdAdjYuTYuAdjYd + 300*traceYuAdjYdTYdAdjYu + 1800*
      traceYuAdjYuTYuAdjYu - 80*traceAdjYuTYu*Sqr(g1) + 80*MassB*traceYuAdjYu*
      Sqr(g1) - 150*traceAdjfuTfu*Sqr(g1p) + 30*traceAdjYuTYu*Sqr(g1p) + 150*
      MassBp*tracefuAdjfu*Sqr(g1p) - 30*MassBp*traceYuAdjYu*Sqr(g1p) + 27*MassB
      *Sqr(g1)*Sqr(g1p) + 27*MassBp*Sqr(g1)*Sqr(g1p) + 180*MassB*Sqr(g1)*Sqr(g2
      ) + 180*MassWB*Sqr(g1)*Sqr(g2) + 195*MassBp*Sqr(g1p)*Sqr(g2) + 195*MassWB
      *Sqr(g1p)*Sqr(g2) - 1600*traceAdjYuTYu*Sqr(g3) + 1600*MassG*traceYuAdjYu*
      Sqr(g3) + 19*Power(g1p,4)*MassBp*Sqr(QS) + 600*Lambdax*Sqr(Conj(Lambdax))
      *TLambdax + 5*Conj(Lambdax)*((20*tracefdAdjfd + 60*traceKappaAdjKappa +
      40*traceLambda12AdjLambda12 + 60*traceYdAdjYd + 20*traceYeAdjYe + 20*
      AbsSqr(Sigmax) - 5*Sqr(g1p) - Sqr(g1p)*Sqr(QS))*TLambdax + Lambdax*(20*
      traceAdjfdTfd + 60*traceAdjKappaTKappa + 40*traceAdjLambda12TLambda12 +
      60*traceAdjYdTYd + 20*traceAdjYeTYe + 5*MassBp*Sqr(g1p) + MassBp*Sqr(g1p)
      *Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) - 0.8*(5*traceAdjfdTfd + 15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 3*MassB*Sqr(g1) - 3*MassBp*Sqr(g1p) +
      15*MassWB*Sqr(g2) + 5*Conj(Lambdax)*TLambdax)*(fd*fd.adjoint()*fu) - 2*
      tracefdAdjfd*(fd*fd.adjoint()*Tfu) - 6*traceYdAdjYd*(fd*fd.adjoint()*Tfu)
      - 2*traceYeAdjYe*(fd*fd.adjoint()*Tfu) - 2*AbsSqr(Lambdax)*(fd*
      fd.adjoint()*Tfu) + 1.2*Sqr(g1)*(fd*fd.adjoint()*Tfu) - 1.2*Sqr(g1p)*(fd*
      fd.adjoint()*Tfu) + 6*Sqr(g2)*(fd*fd.adjoint()*Tfu) - 6*traceAdjfuTfu*(fu
      *fu.adjoint()*fu) - 18*traceAdjYuTYu*(fu*fu.adjoint()*fu) - 2.4*MassB*Sqr
      (g1)*(fu*fu.adjoint()*fu) + 0.4*MassBp*Sqr(g1p)*(fu*fu.adjoint()*fu) - 12
      *MassWB*Sqr(g2)*(fu*fu.adjoint()*fu) - 6*Conj(Lambdax)*TLambdax*(fu*
      fu.adjoint()*fu) - 4*tracefuAdjfu*(fu*fu.adjoint()*Tfu) - 12*traceYuAdjYu
      *(fu*fu.adjoint()*Tfu) - 4*AbsSqr(Lambdax)*(fu*fu.adjoint()*Tfu) + 1.2*
      Sqr(g1)*(fu*fu.adjoint()*Tfu) + 0.8*Sqr(g1p)*(fu*fu.adjoint()*Tfu) + 6*
      Sqr(g2)*(fu*fu.adjoint()*Tfu) - 6*traceAdjgDTgD*(fu*hE.adjoint()*hE) - 2*
      traceAdjhEThE*(fu*hE.adjoint()*hE) - 2.4*MassB*Sqr(g1)*(fu*hE.adjoint()*
      hE) + 0.4*MassBp*Sqr(g1p)*(fu*hE.adjoint()*hE) - 2*Conj(SigmaL)*TSigmaL*(
      fu*hE.adjoint()*hE) - 6*tracegDAdjgD*(fu*hE.adjoint()*ThE) - 2*
      tracehEAdjhE*(fu*hE.adjoint()*ThE) - 2*AbsSqr(SigmaL)*(fu*hE.adjoint()*
      ThE) + 2.4*Sqr(g1)*(fu*hE.adjoint()*ThE) - 0.4*Sqr(g1p)*(fu*hE.adjoint()*
      ThE) - 6*traceAdjKappaTKappa*(fu*(Lambda12).adjoint()*Lambda12) - 4*
      traceAdjLambda12TLambda12*(fu*(Lambda12).adjoint()*Lambda12) + 0.5*MassBp
      *Sqr(g1p)*(fu*(Lambda12).adjoint()*Lambda12) - 0.1*MassBp*Sqr(g1p)*Sqr(QS
      )*(fu*(Lambda12).adjoint()*Lambda12) - 4*Conj(Lambdax)*TLambdax*(fu*(
      Lambda12).adjoint()*Lambda12) - 2*Conj(Sigmax)*TSigmax*(fu*(Lambda12)
      .adjoint()*Lambda12) - 6*traceKappaAdjKappa*(fu*(Lambda12).adjoint()*
      TLambda12) - 4*traceLambda12AdjLambda12*(fu*(Lambda12).adjoint()*
      TLambda12) - 4*AbsSqr(Lambdax)*(fu*(Lambda12).adjoint()*TLambda12) - 2*
      AbsSqr(Sigmax)*(fu*(Lambda12).adjoint()*TLambda12) - 0.5*Sqr(g1p)*(fu*(
      Lambda12).adjoint()*TLambda12) + 0.1*Sqr(g1p)*Sqr(QS)*(fu*(Lambda12)
      .adjoint()*TLambda12) - 4*tracefdAdjfd*(Tfd*fd.adjoint()*fu) - 12*
      traceYdAdjYd*(Tfd*fd.adjoint()*fu) - 4*traceYeAdjYe*(Tfd*fd.adjoint()*fu)
      - 4*AbsSqr(Lambdax)*(Tfd*fd.adjoint()*fu) + 2.4*Sqr(g1)*(Tfd*fd.adjoint(
      )*fu) - 2.4*Sqr(g1p)*(Tfd*fd.adjoint()*fu) + 12*Sqr(g2)*(Tfd*fd.adjoint()
      *fu) - 5*tracefuAdjfu*(Tfu*fu.adjoint()*fu) - 15*traceYuAdjYu*(Tfu*
      fu.adjoint()*fu) - 5*AbsSqr(Lambdax)*(Tfu*fu.adjoint()*fu) + 2.4*Sqr(g1)*
      (Tfu*fu.adjoint()*fu) - 1.4*Sqr(g1p)*(Tfu*fu.adjoint()*fu) + 12*Sqr(g2)*(
      Tfu*fu.adjoint()*fu) - 3*tracegDAdjgD*(Tfu*hE.adjoint()*hE) -
      tracehEAdjhE*(Tfu*hE.adjoint()*hE) - AbsSqr(SigmaL)*(Tfu*hE.adjoint()*hE)
      + 1.2*Sqr(g1)*(Tfu*hE.adjoint()*hE) - 0.2*Sqr(g1p)*(Tfu*hE.adjoint()*hE)
      - 3*traceKappaAdjKappa*(Tfu*(Lambda12).adjoint()*Lambda12) - 2*
      traceLambda12AdjLambda12*(Tfu*(Lambda12).adjoint()*Lambda12) - 2*AbsSqr(
      Lambdax)*(Tfu*(Lambda12).adjoint()*Lambda12) - AbsSqr(Sigmax)*(Tfu*(
      Lambda12).adjoint()*Lambda12) - 0.25*Sqr(g1p)*(Tfu*(Lambda12).adjoint()*
      Lambda12) + 0.05*Sqr(g1p)*Sqr(QS)*(Tfu*(Lambda12).adjoint()*Lambda12) - 2
      *(fd*fd.adjoint()*fd*fd.adjoint()*Tfu) - 4*(fd*fd.adjoint()*Tfd*
      fd.adjoint()*fu) - 2*(fd*Lambda12.conjugate()*(Lambda12).transpose()*
      fd.adjoint()*Tfu) - 4*(fd*Lambda12.conjugate()*(TLambda12).transpose()*
      fd.adjoint()*fu) - 4*(fu*fu.adjoint()*fd*fd.adjoint()*Tfu) - 6*(fu*
      fu.adjoint()*fu*fu.adjoint()*Tfu) - 4*(fu*fu.adjoint()*Tfd*fd.adjoint()*
      fu) - 8*(fu*fu.adjoint()*Tfu*fu.adjoint()*fu) - 2*(fu*hE.adjoint()*hE*
      fu.adjoint()*Tfu) - 4*(fu*hE.adjoint()*hE*hE.adjoint()*ThE) - 4*(fu*
      hE.adjoint()*Ye*Ye.adjoint()*ThE) - 4*(fu*hE.adjoint()*ThE*fu.adjoint()*
      fu) - 4*(fu*hE.adjoint()*ThE*hE.adjoint()*hE) - 4*(fu*hE.adjoint()*TYe*
      Ye.adjoint()*hE) - 2*(fu*(Lambda12).adjoint()*Lambda12*fu.adjoint()*Tfu)
      - 2*(fu*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*TLambda12) - 4
      *(fu*(Lambda12).adjoint()*TLambda12*fu.adjoint()*fu) - 2*(fu*(Lambda12)
      .adjoint()*TLambda12*(Lambda12).adjoint()*Lambda12) - 2*(fu*(Lambda12)
      .adjoint()*fd.transpose()*fd.conjugate()*TLambda12) - 2*(fu*(Lambda12)
      .adjoint()*(Tfd).transpose()*fd.conjugate()*Lambda12) - 4*(Tfd*fd.adjoint
      ()*fd*fd.adjoint()*fu) - 4*(Tfd*Lambda12.conjugate()*(Lambda12).transpose
      ()*fd.adjoint()*fu) - 2*(Tfu*fu.adjoint()*fd*fd.adjoint()*fu) - 6*(Tfu*
      fu.adjoint()*fu*fu.adjoint()*fu) - 4*(Tfu*hE.adjoint()*hE*fu.adjoint()*fu
      ) - 2*(Tfu*hE.adjoint()*hE*hE.adjoint()*hE) - 2*(Tfu*hE.adjoint()*Ye*
      Ye.adjoint()*hE) - 4*(Tfu*(Lambda12).adjoint()*Lambda12*fu.adjoint()*fu)
      - Tfu*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12 - Tfu*(
      Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12)).real();


   return beta_Tfu;
}

} // namespace flexiblesusy
