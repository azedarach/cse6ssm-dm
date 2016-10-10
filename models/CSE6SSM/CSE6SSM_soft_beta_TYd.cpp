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

// File generated at Wed 3 Jun 2015 23:42:58

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(tracefdAdjfd*TYd + 3*traceYdAdjYd*TYd +
      traceYeAdjYe*TYd + AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd -
      0.7*Sqr(g1p)*TYd - 3*Sqr(g2)*TYd - 5.333333333333333*Sqr(g3)*TYd + Yd*(2
      *traceAdjfdTfd + 6*traceAdjYdTYd + 2*traceAdjYeTYe + 0.9333333333333333*
      MassB*Sqr(g1) + 1.4*MassBp*Sqr(g1p) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 2*Conj(Lambdax)*TLambdax) + 4*(Yd*
      Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 2*(Yd*gD.conjugate()*(TgD)
      .transpose()) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu + TYd*
      gD.conjugate()*gD.transpose())).real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double tracefdAdjfdTfdAdjfd = TRACE_STRUCT.tracefdAdjfdTfdAdjfd;
   const double tracefdAdjfdTfuAdjfu = TRACE_STRUCT.tracefdAdjfdTfuAdjfu;
   const double tracefuAdjfuTfdAdjfd = TRACE_STRUCT.tracefuAdjfuTfdAdjfd;
   const double tracehEAdjhETYeAdjYe = TRACE_STRUCT.tracehEAdjhETYeAdjYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeThEAdjhE = TRACE_STRUCT.traceYeAdjYeThEAdjhE;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceAdjfdTfdconjLambda12TpLambda12 =
      TRACE_STRUCT.traceAdjfdTfdconjLambda12TpLambda12;
   const double traceAdjgDTpYdconjYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjYdTgD;
   const double traceAdjYdTYdconjgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjgDTpgD;
   const double traceAdjLambda12TpfdconjfdTLambda12 =
      TRACE_STRUCT.traceAdjLambda12TpfdconjfdTLambda12;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
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


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(4.588888888888889*Power(g1,4)*TYd + 6.825*Power(
      g1p,4)*TYd + 16.5*Power(g2,4)*TYd + 14.222222222222221*Power(g3,4)*TYd -
      3*tracefdAdjfdfdAdjfd*TYd - 2*tracefdAdjfdfuAdjfu*TYd - 3*
      tracegDAdjgDTpYdconjYd*TYd - 2*tracehEAdjhEYeAdjYe*TYd -
      traceLambda12AdjLambda12Tpfdconjfd*TYd - 9*traceYdAdjYdYdAdjYd*TYd - 3*
      traceYdAdjYuYuAdjYd*TYd - 3*traceYeAdjYeYeAdjYe*TYd - tracefuAdjfu*AbsSqr
      (Lambdax)*TYd - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd - 2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYd - 3*traceYuAdjYu*AbsSqr(
      Lambdax)*TYd - AbsSqr(Lambdax)*AbsSqr(Sigmax)*TYd - 0.4*traceYdAdjYd*Sqr(
      g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*TYd + tracefdAdjfd*Sqr(g1p)*TYd - 0.6*
      traceYdAdjYd*Sqr(g1p)*TYd - 0.2*traceYeAdjYe*Sqr(g1p)*TYd - 0.25*AbsSqr(
      Lambdax)*Sqr(g1p)*TYd - 0.23333333333333334*Sqr(g1)*Sqr(g1p)*TYd + Sqr(g1
      )*Sqr(g2)*TYd + 1.5*Sqr(g1p)*Sqr(g2)*TYd + 16*traceYdAdjYd*Sqr(g3)*TYd +
      0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 1.3333333333333333*Sqr(g1p)*Sqr(
      g3)*TYd + 8*Sqr(g2)*Sqr(g3)*TYd + 0.035*Power(g1p,4)*Sqr(QS)*TYd + 0.05*
      AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TYd - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*
      TYd - 0.0022222222222222222*Yd*(8260*Power(g1,4)*MassB + 12285*Power(g1p,
      4)*MassBp + 25600*Power(g3,4)*MassG + 29700*Power(g2,4)*MassWB + 900*
      traceAdjfdTfdconjLambda12TpLambda12 + 2700*traceAdjgDTpYdconjYdTgD + 900*
      traceAdjLambda12TpfdconjfdTLambda12 + 2700*traceAdjYdTYdconjgDTpgD + 5400
      *tracefdAdjfdTfdAdjfd + 1800*tracefdAdjfdTfuAdjfu + 1800*
      tracefuAdjfuTfdAdjfd + 1800*tracehEAdjhETYeAdjYe + 16200*
      traceYdAdjYdTYdAdjYd + 2700*traceYdAdjYuTYuAdjYd + 1800*
      traceYeAdjYeThEAdjhE + 5400*traceYeAdjYeTYeAdjYe + 2700*
      traceYuAdjYdTYdAdjYu + 360*traceAdjYdTYd*Sqr(g1) - 1080*traceAdjYeTYe*Sqr
      (g1) - 360*MassB*traceYdAdjYd*Sqr(g1) + 1080*MassB*traceYeAdjYe*Sqr(g1) -
      900*traceAdjfdTfd*Sqr(g1p) + 540*traceAdjYdTYd*Sqr(g1p) + 180*
      traceAdjYeTYe*Sqr(g1p) + 900*MassBp*tracefdAdjfd*Sqr(g1p) - 540*MassBp*
      traceYdAdjYd*Sqr(g1p) - 180*MassBp*traceYeAdjYe*Sqr(g1p) - 210*MassB*Sqr(
      g1)*Sqr(g1p) - 210*MassBp*Sqr(g1)*Sqr(g1p) + 900*MassB*Sqr(g1)*Sqr(g2) +
      900*MassWB*Sqr(g1)*Sqr(g2) + 1350*MassBp*Sqr(g1p)*Sqr(g2) + 1350*MassWB*
      Sqr(g1p)*Sqr(g2) - 14400*traceAdjYdTYd*Sqr(g3) + 14400*MassG*traceYdAdjYd
      *Sqr(g3) + 800*MassB*Sqr(g1)*Sqr(g3) + 800*MassG*Sqr(g1)*Sqr(g3) + 1200*
      MassBp*Sqr(g1p)*Sqr(g3) + 1200*MassG*Sqr(g1p)*Sqr(g3) + 7200*MassG*Sqr(g2
      )*Sqr(g3) + 7200*MassWB*Sqr(g2)*Sqr(g3) + 63*Power(g1p,4)*MassBp*Sqr(QS)
      + 5400*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 45*Conj(Lambdax)*((20*
      tracefuAdjfu + 60*traceKappaAdjKappa + 40*traceLambda12AdjLambda12 + 60*
      traceYuAdjYu + 20*AbsSqr(Sigmax) + 5*Sqr(g1p) - Sqr(g1p)*Sqr(QS))*
      TLambdax + Lambdax*(20*traceAdjfuTfu + 60*traceAdjKappaTKappa + 40*
      traceAdjLambda12TLambda12 + 60*traceAdjYuTYu - 5*MassBp*Sqr(g1p) + MassBp
      *Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) - 0.4*(15*traceAdjfdTfd +
      45*traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr(g1) + 6*MassBp*Sqr(g1p)
      + 30*MassWB*Sqr(g2) + 15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) -
      4*tracefdAdjfd*(Yd*Yd.adjoint()*TYd) - 12*traceYdAdjYd*(Yd*Yd.adjoint()*
      TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*TYd) - 4*AbsSqr(Lambdax)*(Yd*
      Yd.adjoint()*TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()*TYd) + 1.8*Sqr(g1p)*(Yd*
      Yd.adjoint()*TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd) - 2*traceAdjfuTfu*(Yd
      *Yu.adjoint()*Yu) - 6*traceAdjYuTYu*(Yd*Yu.adjoint()*Yu) - 1.6*MassB*Sqr(
      g1)*(Yd*Yu.adjoint()*Yu) - 0.4*MassBp*Sqr(g1p)*(Yd*Yu.adjoint()*Yu) - 2*
      Conj(Lambdax)*TLambdax*(Yd*Yu.adjoint()*Yu) - 2*tracefuAdjfu*(Yd*
      Yu.adjoint()*TYu) - 6*traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*AbsSqr(
      Lambdax)*(Yd*Yu.adjoint()*TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) + 0.4*
      Sqr(g1p)*(Yd*Yu.adjoint()*TYu) - 6*traceAdjgDTgD*(Yd*gD.conjugate()*
      gD.transpose()) - 2*traceAdjhEThE*(Yd*gD.conjugate()*gD.transpose()) -
      0.8*MassB*Sqr(g1)*(Yd*gD.conjugate()*gD.transpose()) - 1.2*MassBp*Sqr(g1p
      )*(Yd*gD.conjugate()*gD.transpose()) - 2*Conj(SigmaL)*TSigmaL*(Yd*
      gD.conjugate()*gD.transpose()) - 6*tracegDAdjgD*(Yd*gD.conjugate()*(TgD)
      .transpose()) - 2*tracehEAdjhE*(Yd*gD.conjugate()*(TgD).transpose()) - 2*
      AbsSqr(SigmaL)*(Yd*gD.conjugate()*(TgD).transpose()) + 0.8*Sqr(g1)*(Yd*
      gD.conjugate()*(TgD).transpose()) + 1.2*Sqr(g1p)*(Yd*gD.conjugate()*(TgD)
      .transpose()) - 5*tracefdAdjfd*(TYd*Yd.adjoint()*Yd) - 15*traceYdAdjYd*(
      TYd*Yd.adjoint()*Yd) - 5*traceYeAdjYe*(TYd*Yd.adjoint()*Yd) - 5*AbsSqr(
      Lambdax)*(TYd*Yd.adjoint()*Yd) + 1.2*Sqr(g1)*(TYd*Yd.adjoint()*Yd) + 1.8*
      Sqr(g1p)*(TYd*Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*Yd) -
      tracefuAdjfu*(TYd*Yu.adjoint()*Yu) - 3*traceYuAdjYu*(TYd*Yu.adjoint()*Yu)
      - AbsSqr(Lambdax)*(TYd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(TYd*Yu.adjoint()*
      Yu) + 0.2*Sqr(g1p)*(TYd*Yu.adjoint()*Yu) - 3*tracegDAdjgD*(TYd*
      gD.conjugate()*gD.transpose()) - tracehEAdjhE*(TYd*gD.conjugate()*
      gD.transpose()) - AbsSqr(SigmaL)*(TYd*gD.conjugate()*gD.transpose()) +
      0.4*Sqr(g1)*(TYd*gD.conjugate()*gD.transpose()) + 0.6*Sqr(g1p)*(TYd*
      gD.conjugate()*gD.transpose()) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd)
      - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*
      Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 2*(Yd*gD.conjugate()*gD.transpose()*Yd.adjoint()*TYd) - 4*(Yd*
      gD.conjugate()*gD.transpose()*gD.conjugate()*(TgD).transpose()) - 2*(Yd*
      gD.conjugate()*(Kappa).transpose()*Kappa.conjugate()*(TgD).transpose()) -
      4*(Yd*gD.conjugate()*(TgD).transpose()*Yd.adjoint()*Yd) - 4*(Yd*
      gD.conjugate()*(TgD).transpose()*gD.conjugate()*gD.transpose()) - 2*(Yd*
      gD.conjugate()*(TKappa).transpose()*Kappa.conjugate()*gD.transpose()) - 6
      *(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 4*(TYd*
      gD.conjugate()*gD.transpose()*Yd.adjoint()*Yd) - 2*(TYd*gD.conjugate()*
      gD.transpose()*gD.conjugate()*gD.transpose()) - TYd*gD.conjugate()*(Kappa
      ).transpose()*Kappa.conjugate()*gD.transpose())).real();


   return beta_TYd;
}

} // namespace flexiblesusy
