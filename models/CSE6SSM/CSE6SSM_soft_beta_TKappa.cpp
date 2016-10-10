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

// File generated at Wed 3 Jun 2015 23:43:06

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TKappa_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (oneOver16PiSqr*(3*traceKappaAdjKappa*TKappa + 2*
      traceLambda12AdjLambda12*TKappa + 2*AbsSqr(Lambdax)*TKappa + AbsSqr(
      Sigmax)*TKappa - 0.26666666666666666*Sqr(g1)*TKappa - 0.65*Sqr(g1p)*
      TKappa - 5.333333333333333*Sqr(g3)*TKappa - 0.05*Sqr(g1p)*Sqr(QS)*TKappa
      + Kappa*(6*traceAdjKappaTKappa + 4*traceAdjLambda12TLambda12 +
      0.5333333333333333*MassB*Sqr(g1) + 1.3*MassBp*Sqr(g1p) +
      10.666666666666666*MassG*Sqr(g3) + 0.1*MassBp*Sqr(g1p)*Sqr(QS) + 4*Conj(
      Lambdax)*TLambdax + 2*Conj(Sigmax)*TSigmax) + 4*(Kappa*gD.adjoint()*TgD)
      + 3*(Kappa*(Kappa).adjoint()*TKappa) + 2*(TKappa*gD.adjoint()*gD) + 3*(
      TKappa*(Kappa).adjoint()*Kappa))).real();


   return beta_TKappa;
}

/**
 * Calculates the two-loop beta function of TKappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TKappa_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (twoLoop*(2.5955555555555554*Power(g1,4)*TKappa + 6.3525
      *Power(g1p,4)*TKappa + 14.222222222222221*Power(g3,4)*TKappa + 0.005*
      Power(g1p,4)*Power(QS,4)*TKappa - 2*tracefuAdjLambda12Lambda12Adjfu*
      TKappa - 6*tracegDAdjKappaKappaAdjgD*TKappa - 2*
      tracehEAdjLambda12Lambda12AdjhE*TKappa - 6*
      traceKappaAdjKappaKappaAdjKappa*TKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TKappa - 2*
      traceLambda12AdjLambda12Tpfdconjfd*TKappa - 2*tracefdAdjfd*AbsSqr(Lambdax
      )*TKappa - 2*tracefuAdjfu*AbsSqr(Lambdax)*TKappa - 6*traceYdAdjYd*AbsSqr(
      Lambdax)*TKappa - 2*traceYeAdjYe*AbsSqr(Lambdax)*TKappa - 6*traceYuAdjYu*
      AbsSqr(Lambdax)*TKappa - 2*AbsSqr(KappaPr)*AbsSqr(Sigmax)*TKappa - 2*
      AbsSqr(Sigmax)*AbsSqr(SigmaL)*TKappa + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TKappa + 1.2*traceLambda12AdjLambda12*Sqr(g1)*TKappa + 1.2*AbsSqr(Lambdax
      )*Sqr(g1)*TKappa + 1.95*traceKappaAdjKappa*Sqr(g1p)*TKappa + 1.3*
      traceLambda12AdjLambda12*Sqr(g1p)*TKappa + 1.3*AbsSqr(Lambdax)*Sqr(g1p)*
      TKappa + 0.25333333333333335*Sqr(g1)*Sqr(g1p)*TKappa + 6*
      traceLambda12AdjLambda12*Sqr(g2)*TKappa + 6*AbsSqr(Lambdax)*Sqr(g2)*
      TKappa + 16*traceKappaAdjKappa*Sqr(g3)*TKappa + 1.4222222222222223*Sqr(g1
      )*Sqr(g3)*TKappa + 3.466666666666667*Sqr(g1p)*Sqr(g3)*TKappa + 0.5025*
      Power(g1p,4)*Sqr(QS)*TKappa - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*
      TKappa - 0.1*traceLambda12AdjLambda12*Sqr(g1p)*Sqr(QS)*TKappa - 0.1*
      AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TKappa - 4*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TKappa - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)*TKappa -
      0.0011111111111111111*Kappa*(9344*Power(g1,4)*MassB + 22869*Power(g1p,4)*
      MassBp + 51200*Power(g3,4)*MassG + 18*Power(g1p,4)*MassBp*Power(QS,4) +
      3600*traceAdjfdTfdconjLambda12TpLambda12 + 3600*
      traceAdjLambda12TpfdconjfdTLambda12 + 3600*
      tracefuAdjLambda12TLambda12Adjfu + 10800*tracegDAdjKappaTKappaAdjgD +
      3600*tracehEAdjLambda12TLambda12AdjhE + 10800*traceKappaAdjgDTgDAdjKappa
      + 21600*traceKappaAdjKappaTKappaAdjKappa + 3600*
      traceLambda12AdjfuTfuAdjLambda12 + 3600*traceLambda12AdjhEThEAdjLambda12
      + 14400*traceLambda12AdjLambda12TLambda12AdjLambda12 - 1440*
      traceAdjKappaTKappa*Sqr(g1) - 2160*traceAdjLambda12TLambda12*Sqr(g1) +
      1440*MassB*traceKappaAdjKappa*Sqr(g1) + 2160*MassB*
      traceLambda12AdjLambda12*Sqr(g1) - 3510*traceAdjKappaTKappa*Sqr(g1p) -
      2340*traceAdjLambda12TLambda12*Sqr(g1p) + 3510*MassBp*traceKappaAdjKappa*
      Sqr(g1p) + 2340*MassBp*traceLambda12AdjLambda12*Sqr(g1p) + 456*MassB*Sqr(
      g1)*Sqr(g1p) + 456*MassBp*Sqr(g1)*Sqr(g1p) - 10800*
      traceAdjLambda12TLambda12*Sqr(g2) + 10800*MassWB*traceLambda12AdjLambda12
      *Sqr(g2) - 28800*traceAdjKappaTKappa*Sqr(g3) + 28800*MassG*
      traceKappaAdjKappa*Sqr(g3) + 2560*MassB*Sqr(g1)*Sqr(g3) + 2560*MassG*Sqr(
      g1)*Sqr(g3) + 6240*MassBp*Sqr(g1p)*Sqr(g3) + 6240*MassG*Sqr(g1p)*Sqr(g3)
      + 1809*Power(g1p,4)*MassBp*Sqr(QS) + 270*traceAdjKappaTKappa*Sqr(g1p)*Sqr
      (QS) + 180*traceAdjLambda12TLambda12*Sqr(g1p)*Sqr(QS) - 270*MassBp*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) - 180*MassBp*traceLambda12AdjLambda12
      *Sqr(g1p)*Sqr(QS) + 14400*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 180*Conj(
      Lambdax)*(Lambdax*(20*traceAdjfdTfd + 20*traceAdjfuTfu + 60*traceAdjYdTYd
      + 20*traceAdjYeTYe + 60*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr
      (g1p) + 60*MassWB*Sqr(g2) - MassBp*Sqr(g1p)*Sqr(QS)) + (20*tracefdAdjfd +
      20*tracefuAdjfu + 60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu -
      12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) + Sqr(g1p)*Sqr(QS))*TLambdax) +
      3600*AbsSqr(SigmaL)*Conj(Sigmax)*TSigmax + 7200*Sigmax*Sqr(Conj(Sigmax))*
      TSigmax + 3600*Conj(KappaPr)*Conj(Sigmax)*(Sigmax*TKappaPr + KappaPr*
      TSigmax) + 3600*AbsSqr(Sigmax)*Conj(SigmaL)*TSigmaL) - 0.8*(15*
      traceAdjgDTgD + 5*traceAdjhEThE + MassB*Sqr(g1) - MassBp*Sqr(g1p) + 15*
      MassWB*Sqr(g2) + 5*Conj(SigmaL)*TSigmaL)*(Kappa*gD.adjoint()*gD) - 12*
      tracegDAdjgD*(Kappa*gD.adjoint()*TgD) - 4*tracehEAdjhE*(Kappa*gD.adjoint(
      )*TgD) - 4*AbsSqr(SigmaL)*(Kappa*gD.adjoint()*TgD) + 0.8*Sqr(g1)*(Kappa*
      gD.adjoint()*TgD) - 0.8*Sqr(g1p)*(Kappa*gD.adjoint()*TgD) + 12*Sqr(g2)*(
      Kappa*gD.adjoint()*TgD) - 12*traceAdjKappaTKappa*(Kappa*(Kappa).adjoint()
      *Kappa) - 8*traceAdjLambda12TLambda12*(Kappa*(Kappa).adjoint()*Kappa) -
      0.2*MassBp*Sqr(g1p)*Sqr(QS)*(Kappa*(Kappa).adjoint()*Kappa) - 8*Conj(
      Lambdax)*TLambdax*(Kappa*(Kappa).adjoint()*Kappa) - 4*Conj(Sigmax)*
      TSigmax*(Kappa*(Kappa).adjoint()*Kappa) - 9*traceKappaAdjKappa*(Kappa*(
      Kappa).adjoint()*TKappa) - 6*traceLambda12AdjLambda12*(Kappa*(Kappa)
      .adjoint()*TKappa) - 6*AbsSqr(Lambdax)*(Kappa*(Kappa).adjoint()*TKappa) -
      3*AbsSqr(Sigmax)*(Kappa*(Kappa).adjoint()*TKappa) - 0.25*Sqr(g1p)*(Kappa
      *(Kappa).adjoint()*TKappa) + 0.15*Sqr(g1p)*Sqr(QS)*(Kappa*(Kappa).adjoint
      ()*TKappa) - 6*tracegDAdjgD*(TKappa*gD.adjoint()*gD) - 2*tracehEAdjhE*(
      TKappa*gD.adjoint()*gD) - 2*AbsSqr(SigmaL)*(TKappa*gD.adjoint()*gD) + 0.4
      *Sqr(g1)*(TKappa*gD.adjoint()*gD) - 0.4*Sqr(g1p)*(TKappa*gD.adjoint()*gD)
      + 6*Sqr(g2)*(TKappa*gD.adjoint()*gD) - 9*traceKappaAdjKappa*(TKappa*(
      Kappa).adjoint()*Kappa) - 6*traceLambda12AdjLambda12*(TKappa*(Kappa)
      .adjoint()*Kappa) - 6*AbsSqr(Lambdax)*(TKappa*(Kappa).adjoint()*Kappa) -
      3*AbsSqr(Sigmax)*(TKappa*(Kappa).adjoint()*Kappa) + 0.25*Sqr(g1p)*(TKappa
      *(Kappa).adjoint()*Kappa) + 0.15*Sqr(g1p)*Sqr(QS)*(TKappa*(Kappa).adjoint
      ()*Kappa) - 4*(Kappa*gD.adjoint()*gD*gD.adjoint()*TgD) - 2*(Kappa*
      gD.adjoint()*gD*(Kappa).adjoint()*TKappa) - 4*(Kappa*gD.adjoint()*TgD*
      gD.adjoint()*gD) - 4*(Kappa*gD.adjoint()*TgD*(Kappa).adjoint()*Kappa) - 4
      *(Kappa*gD.adjoint()*Yd.transpose()*Yd.conjugate()*TgD) - 4*(Kappa*
      gD.adjoint()*Yu.transpose()*Yu.conjugate()*TgD) - 4*(Kappa*gD.adjoint()*(
      TYd).transpose()*Yd.conjugate()*gD) - 4*(Kappa*gD.adjoint()*(TYu)
      .transpose()*Yu.conjugate()*gD) - 3*(Kappa*(Kappa).adjoint()*Kappa*(Kappa
      ).adjoint()*TKappa) - 4*(Kappa*(Kappa).adjoint()*TKappa*(Kappa).adjoint()
      *Kappa) - 2*(TKappa*gD.adjoint()*gD*gD.adjoint()*gD) - 4*(TKappa*
      gD.adjoint()*gD*(Kappa).adjoint()*Kappa) - 2*(TKappa*gD.adjoint()*
      Yd.transpose()*Yd.conjugate()*gD) - 2*(TKappa*gD.adjoint()*Yu.transpose()
      *Yu.conjugate()*gD) - 3*(TKappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()
      *Kappa))).real();


   return beta_TKappa;
}

} // namespace flexiblesusy
