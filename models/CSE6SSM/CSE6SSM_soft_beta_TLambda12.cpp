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

// File generated at Wed 3 Jun 2015 23:43:08

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_TLambda12_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (oneOver16PiSqr*(3*traceKappaAdjKappa*TLambda12 + 2*
      traceLambda12AdjLambda12*TLambda12 + 2*AbsSqr(Lambdax)*TLambda12 + AbsSqr
      (Sigmax)*TLambda12 - 0.6*Sqr(g1)*TLambda12 - 0.65*Sqr(g1p)*TLambda12 - 3*
      Sqr(g2)*TLambda12 - 0.05*Sqr(g1p)*Sqr(QS)*TLambda12 + Lambda12*(6*
      traceAdjKappaTKappa + 4*traceAdjLambda12TLambda12 + 1.2*MassB*Sqr(g1) +
      1.3*MassBp*Sqr(g1p) + 6*MassWB*Sqr(g2) + 0.1*MassBp*Sqr(g1p)*Sqr(QS) + 4*
      Conj(Lambdax)*TLambdax + 2*Conj(Sigmax)*TSigmax) + 2*(Lambda12*fu.adjoint
      ()*Tfu) + 2*(Lambda12*hE.adjoint()*ThE) + 3*(Lambda12*(Lambda12).adjoint(
      )*TLambda12) + TLambda12*fu.adjoint()*fu + TLambda12*hE.adjoint()*hE + 3*
      (TLambda12*(Lambda12).adjoint()*Lambda12) + fd.transpose()*fd.conjugate()
      *TLambda12 + 2*((Tfd).transpose()*fd.conjugate()*Lambda12))).real();


   return beta_TLambda12;
}

/**
 * Calculates the two-loop beta function of TLambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_TLambda12_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double tracefuAdjLambda12TLambda12Adjfu =
      TRACE_STRUCT.tracefuAdjLambda12TLambda12Adjfu;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE;
   const double traceKappaAdjgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjKappa;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjfuTfuAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjfuTfuAdjLambda12;
   const double traceLambda12AdjhEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceAdjfdTfdconjLambda12TpLambda12 =
      TRACE_STRUCT.traceAdjfdTfdconjLambda12TpLambda12;
   const double traceAdjLambda12TpfdconjfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjfdTLambda12;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
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


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (twoLoop*(5.94*Power(g1,4)*TLambda12 + 6.3525*Power(
      g1p,4)*TLambda12 + 16.5*Power(g2,4)*TLambda12 + 0.005*Power(g1p,4)*Power(
      QS,4)*TLambda12 - 2*tracefuAdjLambda12Lambda12Adjfu*TLambda12 - 6*
      tracegDAdjKappaKappaAdjgD*TLambda12 - 2*tracehEAdjLambda12Lambda12AdjhE*
      TLambda12 - 6*traceKappaAdjKappaKappaAdjKappa*TLambda12 - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TLambda12 - 2*
      traceLambda12AdjLambda12Tpfdconjfd*TLambda12 - 2*tracefdAdjfd*AbsSqr(
      Lambdax)*TLambda12 - 2*tracefuAdjfu*AbsSqr(Lambdax)*TLambda12 - 6*
      traceYdAdjYd*AbsSqr(Lambdax)*TLambda12 - 2*traceYeAdjYe*AbsSqr(Lambdax)*
      TLambda12 - 6*traceYuAdjYu*AbsSqr(Lambdax)*TLambda12 - 2*AbsSqr(KappaPr)*
      AbsSqr(Sigmax)*TLambda12 - 2*AbsSqr(Sigmax)*AbsSqr(SigmaL)*TLambda12 +
      0.8*traceKappaAdjKappa*Sqr(g1)*TLambda12 + 1.2*traceLambda12AdjLambda12*
      Sqr(g1)*TLambda12 + 1.2*AbsSqr(Lambdax)*Sqr(g1)*TLambda12 + 1.95*
      traceKappaAdjKappa*Sqr(g1p)*TLambda12 + 1.3*traceLambda12AdjLambda12*Sqr(
      g1p)*TLambda12 + 1.3*AbsSqr(Lambdax)*Sqr(g1p)*TLambda12 + 0.27*Sqr(g1)*
      Sqr(g1p)*TLambda12 + 6*traceLambda12AdjLambda12*Sqr(g2)*TLambda12 + 6*
      AbsSqr(Lambdax)*Sqr(g2)*TLambda12 + 1.8*Sqr(g1)*Sqr(g2)*TLambda12 + 1.95*
      Sqr(g1p)*Sqr(g2)*TLambda12 + 16*traceKappaAdjKappa*Sqr(g3)*TLambda12 +
      0.5025*Power(g1p,4)*Sqr(QS)*TLambda12 - 0.15*traceKappaAdjKappa*Sqr(g1p)*
      Sqr(QS)*TLambda12 - 0.1*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS)*
      TLambda12 - 0.1*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TLambda12 - 4*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TLambda12 - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)*
      TLambda12 - 0.01*Lambda12*(2376*Power(g1,4)*MassB + 2541*Power(g1p,4)*
      MassBp + 6600*Power(g2,4)*MassWB + 2*Power(g1p,4)*MassBp*Power(QS,4) +
      400*traceAdjfdTfdconjLambda12TpLambda12 + 400*
      traceAdjLambda12TpfdconjfdTLambda12 + 400*
      tracefuAdjLambda12TLambda12Adjfu + 1200*tracegDAdjKappaTKappaAdjgD + 400*
      tracehEAdjLambda12TLambda12AdjhE + 1200*traceKappaAdjgDTgDAdjKappa + 2400
      *traceKappaAdjKappaTKappaAdjKappa + 400*traceLambda12AdjfuTfuAdjLambda12
      + 400*traceLambda12AdjhEThEAdjLambda12 + 1600*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 160*traceAdjKappaTKappa*
      Sqr(g1) - 240*traceAdjLambda12TLambda12*Sqr(g1) + 160*MassB*
      traceKappaAdjKappa*Sqr(g1) + 240*MassB*traceLambda12AdjLambda12*Sqr(g1) -
      390*traceAdjKappaTKappa*Sqr(g1p) - 260*traceAdjLambda12TLambda12*Sqr(g1p
      ) + 390*MassBp*traceKappaAdjKappa*Sqr(g1p) + 260*MassBp*
      traceLambda12AdjLambda12*Sqr(g1p) + 54*MassB*Sqr(g1)*Sqr(g1p) + 54*MassBp
      *Sqr(g1)*Sqr(g1p) - 1200*traceAdjLambda12TLambda12*Sqr(g2) + 1200*MassWB*
      traceLambda12AdjLambda12*Sqr(g2) + 360*MassB*Sqr(g1)*Sqr(g2) + 360*MassWB
      *Sqr(g1)*Sqr(g2) + 390*MassBp*Sqr(g1p)*Sqr(g2) + 390*MassWB*Sqr(g1p)*Sqr(
      g2) - 3200*traceAdjKappaTKappa*Sqr(g3) + 3200*MassG*traceKappaAdjKappa*
      Sqr(g3) + 201*Power(g1p,4)*MassBp*Sqr(QS) + 30*traceAdjKappaTKappa*Sqr(
      g1p)*Sqr(QS) + 20*traceAdjLambda12TLambda12*Sqr(g1p)*Sqr(QS) - 30*MassBp*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) - 20*MassBp*traceLambda12AdjLambda12*
      Sqr(g1p)*Sqr(QS) + 1600*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 20*Conj(
      Lambdax)*(Lambdax*(20*traceAdjfdTfd + 20*traceAdjfuTfu + 60*traceAdjYdTYd
      + 20*traceAdjYeTYe + 60*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr
      (g1p) + 60*MassWB*Sqr(g2) - MassBp*Sqr(g1p)*Sqr(QS)) + (20*tracefdAdjfd +
      20*tracefuAdjfu + 60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu -
      12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) + Sqr(g1p)*Sqr(QS))*TLambdax) + 400
      *AbsSqr(SigmaL)*Conj(Sigmax)*TSigmax + 800*Sigmax*Sqr(Conj(Sigmax))*
      TSigmax + 400*Conj(KappaPr)*Conj(Sigmax)*(Sigmax*TKappaPr + KappaPr*
      TSigmax) + 400*AbsSqr(Sigmax)*Conj(SigmaL)*TSigmaL) - 2*(traceAdjfuTfu +
      3*traceAdjYuTYu + MassBp*Sqr(g1p) + Conj(Lambdax)*TLambdax)*(Lambda12*
      fu.adjoint()*fu) - 2*tracefuAdjfu*(Lambda12*fu.adjoint()*Tfu) - 6*
      traceYuAdjYu*(Lambda12*fu.adjoint()*Tfu) - 2*AbsSqr(Lambdax)*(Lambda12*
      fu.adjoint()*Tfu) + 2*Sqr(g1p)*(Lambda12*fu.adjoint()*Tfu) - 6*
      traceAdjgDTgD*(Lambda12*hE.adjoint()*hE) - 2*traceAdjhEThE*(Lambda12*
      hE.adjoint()*hE) - 2.4*MassB*Sqr(g1)*(Lambda12*hE.adjoint()*hE) + 0.4*
      MassBp*Sqr(g1p)*(Lambda12*hE.adjoint()*hE) - 2*Conj(SigmaL)*TSigmaL*(
      Lambda12*hE.adjoint()*hE) - 6*tracegDAdjgD*(Lambda12*hE.adjoint()*ThE) -
      2*tracehEAdjhE*(Lambda12*hE.adjoint()*ThE) - 2*AbsSqr(SigmaL)*(Lambda12*
      hE.adjoint()*ThE) + 2.4*Sqr(g1)*(Lambda12*hE.adjoint()*ThE) - 0.4*Sqr(g1p
      )*(Lambda12*hE.adjoint()*ThE) - 12*traceAdjKappaTKappa*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 8*traceAdjLambda12TLambda12*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 0.2*MassBp*Sqr(g1p)*Sqr(QS)*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 8*Conj(Lambdax)*TLambdax*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 4*Conj(Sigmax)*TSigmax*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 9*traceKappaAdjKappa*(Lambda12*(Lambda12)
      .adjoint()*TLambda12) - 6*traceLambda12AdjLambda12*(Lambda12*(Lambda12)
      .adjoint()*TLambda12) - 6*AbsSqr(Lambdax)*(Lambda12*(Lambda12).adjoint()*
      TLambda12) - 3*AbsSqr(Sigmax)*(Lambda12*(Lambda12).adjoint()*TLambda12) -
      0.25*Sqr(g1p)*(Lambda12*(Lambda12).adjoint()*TLambda12) + 0.15*Sqr(g1p)*
      Sqr(QS)*(Lambda12*(Lambda12).adjoint()*TLambda12) - tracefuAdjfu*(
      TLambda12*fu.adjoint()*fu) - 3*traceYuAdjYu*(TLambda12*fu.adjoint()*fu) -
      AbsSqr(Lambdax)*(TLambda12*fu.adjoint()*fu) + Sqr(g1p)*(TLambda12*
      fu.adjoint()*fu) - 3*tracegDAdjgD*(TLambda12*hE.adjoint()*hE) -
      tracehEAdjhE*(TLambda12*hE.adjoint()*hE) - AbsSqr(SigmaL)*(TLambda12*
      hE.adjoint()*hE) + 1.2*Sqr(g1)*(TLambda12*hE.adjoint()*hE) - 0.2*Sqr(g1p)
      *(TLambda12*hE.adjoint()*hE) - 9*traceKappaAdjKappa*(TLambda12*(Lambda12)
      .adjoint()*Lambda12) - 6*traceLambda12AdjLambda12*(TLambda12*(Lambda12)
      .adjoint()*Lambda12) - 6*AbsSqr(Lambdax)*(TLambda12*(Lambda12).adjoint()*
      Lambda12) - 3*AbsSqr(Sigmax)*(TLambda12*(Lambda12).adjoint()*Lambda12) +
      0.25*Sqr(g1p)*(TLambda12*(Lambda12).adjoint()*Lambda12) + 0.15*Sqr(g1p)*
      Sqr(QS)*(TLambda12*(Lambda12).adjoint()*Lambda12) - 2*traceAdjfdTfd*(
      fd.transpose()*fd.conjugate()*Lambda12) - 6*traceAdjYdTYd*(fd.transpose()
      *fd.conjugate()*Lambda12) - 2*traceAdjYeTYe*(fd.transpose()*fd.conjugate(
      )*Lambda12) - 3*MassBp*Sqr(g1p)*(fd.transpose()*fd.conjugate()*Lambda12)
      - 2*Conj(Lambdax)*TLambdax*(fd.transpose()*fd.conjugate()*Lambda12) -
      tracefdAdjfd*(fd.transpose()*fd.conjugate()*TLambda12) - 3*traceYdAdjYd*(
      fd.transpose()*fd.conjugate()*TLambda12) - traceYeAdjYe*(fd.transpose()*
      fd.conjugate()*TLambda12) - AbsSqr(Lambdax)*(fd.transpose()*fd.conjugate(
      )*TLambda12) + 1.5*Sqr(g1p)*(fd.transpose()*fd.conjugate()*TLambda12) - 2
      *tracefdAdjfd*((Tfd).transpose()*fd.conjugate()*Lambda12) - 6*
      traceYdAdjYd*((Tfd).transpose()*fd.conjugate()*Lambda12) - 2*traceYeAdjYe
      *((Tfd).transpose()*fd.conjugate()*Lambda12) - 2*AbsSqr(Lambdax)*((Tfd)
      .transpose()*fd.conjugate()*Lambda12) + 3*Sqr(g1p)*((Tfd).transpose()*
      fd.conjugate()*Lambda12) - 4*(Lambda12*fu.adjoint()*fd*fd.adjoint()*Tfu)
      - 4*(Lambda12*fu.adjoint()*fu*fu.adjoint()*Tfu) - Lambda12*fu.adjoint()*
      fu*(Lambda12).adjoint()*TLambda12 - 4*(Lambda12*fu.adjoint()*Tfd*
      fd.adjoint()*fu) - 4*(Lambda12*fu.adjoint()*Tfu*fu.adjoint()*fu) - 2*(
      Lambda12*fu.adjoint()*Tfu*(Lambda12).adjoint()*Lambda12) - 4*(Lambda12*
      hE.adjoint()*hE*hE.adjoint()*ThE) - Lambda12*hE.adjoint()*hE*(Lambda12)
      .adjoint()*TLambda12 - 4*(Lambda12*hE.adjoint()*Ye*Ye.adjoint()*ThE) - 4*
      (Lambda12*hE.adjoint()*ThE*hE.adjoint()*hE) - 2*(Lambda12*hE.adjoint()*
      ThE*(Lambda12).adjoint()*Lambda12) - 4*(Lambda12*hE.adjoint()*TYe*
      Ye.adjoint()*hE) - 3*(Lambda12*(Lambda12).adjoint()*Lambda12*(Lambda12)
      .adjoint()*TLambda12) - 4*(Lambda12*(Lambda12).adjoint()*TLambda12*(
      Lambda12).adjoint()*Lambda12) - 2*(Lambda12*(Lambda12).adjoint()*
      fd.transpose()*fd.conjugate()*TLambda12) - 2*(Lambda12*(Lambda12).adjoint
      ()*(Tfd).transpose()*fd.conjugate()*Lambda12) - 2*(TLambda12*fu.adjoint()
      *fd*fd.adjoint()*fu) - 2*(TLambda12*fu.adjoint()*fu*fu.adjoint()*fu) - 2*
      (TLambda12*fu.adjoint()*fu*(Lambda12).adjoint()*Lambda12) - 2*(TLambda12*
      hE.adjoint()*hE*hE.adjoint()*hE) - 2*(TLambda12*hE.adjoint()*hE*(Lambda12
      ).adjoint()*Lambda12) - 2*(TLambda12*hE.adjoint()*Ye*Ye.adjoint()*hE) - 3
      *(TLambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12)
      - TLambda12*(Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12 -
      2*(fd.transpose()*fd.conjugate()*fd.transpose()*fd.conjugate()*TLambda12
      ) - 4*(fd.transpose()*fd.conjugate()*(Tfd).transpose()*fd.conjugate()*
      Lambda12) - 2*(fd.transpose()*fu.conjugate()*fu.transpose()*fd.conjugate(
      )*TLambda12) - 4*(fd.transpose()*fu.conjugate()*(Tfu).transpose()*
      fd.conjugate()*Lambda12) - 4*((Tfd).transpose()*fd.conjugate()*
      fd.transpose()*fd.conjugate()*Lambda12) - 4*((Tfd).transpose()*
      fu.conjugate()*fu.transpose()*fd.conjugate()*Lambda12))).real();


   return beta_TLambda12;
}

} // namespace flexiblesusy
