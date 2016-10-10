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

// File generated at Wed 3 Jun 2015 23:43:14

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYu_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(tracefuAdjfu*TYu + 3*traceYuAdjYu*TYu +
      AbsSqr(Lambdax)*TYu - 0.8666666666666667*Sqr(g1)*TYu - 0.3*Sqr(g1p)*TYu -
      3*Sqr(g2)*TYu - 5.333333333333333*Sqr(g3)*TYu + Yu*(2*traceAdjfuTfu + 6*
      traceAdjYuTYu + 1.7333333333333334*MassB*Sqr(g1) + 0.6*MassBp*Sqr(g1p) +
      6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*Conj(Lambdax)*
      TLambdax) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) + 2*(Yu*
      gD.conjugate()*(TgD).transpose()) + TYu*Yd.adjoint()*Yd + 5*(TYu*
      Yu.adjoint()*Yu) + TYu*gD.conjugate()*gD.transpose())).real();


   return beta_TYu;
}

/**
 * Calculates the two-loop beta function of TYu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYu_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(8.695555555555556*Power(g1,4)*TYu + 2.865*Power(
      g1p,4)*TYu + 16.5*Power(g2,4)*TYu + 14.222222222222221*Power(g3,4)*TYu -
      2*tracefdAdjfdfuAdjfu*TYu - 3*tracefuAdjfufuAdjfu*TYu -
      tracefuAdjhEhEAdjfu*TYu - tracefuAdjLambda12Lambda12Adjfu*TYu - 3*
      tracegDAdjgDTpYuconjYu*TYu - 3*traceYdAdjYuYuAdjYd*TYu - 9*
      traceYuAdjYuYuAdjYu*TYu - tracefdAdjfd*AbsSqr(Lambdax)*TYu - 3*
      traceKappaAdjKappa*AbsSqr(Lambdax)*TYu - 2*traceLambda12AdjLambda12*
      AbsSqr(Lambdax)*TYu - 3*traceYdAdjYd*AbsSqr(Lambdax)*TYu - traceYeAdjYe*
      AbsSqr(Lambdax)*TYu - AbsSqr(Lambdax)*AbsSqr(Sigmax)*TYu + 0.8*
      traceYuAdjYu*Sqr(g1)*TYu + 1.5*tracefuAdjfu*Sqr(g1p)*TYu - 0.3*
      traceYuAdjYu*Sqr(g1p)*TYu + 0.25*AbsSqr(Lambdax)*Sqr(g1p)*TYu +
      0.5366666666666666*Sqr(g1)*Sqr(g1p)*TYu + Sqr(g1)*Sqr(g2)*TYu + 0.75*Sqr(
      g1p)*Sqr(g2)*TYu + 16*traceYuAdjYu*Sqr(g3)*TYu + 3.022222222222222*Sqr(g1
      )*Sqr(g3)*TYu + 0.5333333333333333*Sqr(g1p)*Sqr(g3)*TYu + 8*Sqr(g2)*Sqr(
      g3)*TYu + 0.015*Power(g1p,4)*Sqr(QS)*TYu + 0.05*AbsSqr(Lambdax)*Sqr(g1p)*
      Sqr(QS)*TYu - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYu -
      0.0022222222222222222*Yu*(15652*Power(g1,4)*MassB + 5157*Power(g1p,4)*
      MassBp + 25600*Power(g3,4)*MassG + 29700*Power(g2,4)*MassWB + 2700*
      traceAdjgDTpYuconjYuTgD + 2700*traceAdjYuTYuconjgDTpgD + 1800*
      tracefdAdjfdTfuAdjfu + 1800*tracefuAdjfuTfdAdjfd + 5400*
      tracefuAdjfuTfuAdjfu + 900*tracefuAdjhEThEAdjfu + 900*
      tracefuAdjLambda12TLambda12Adjfu + 900*tracehEAdjfuTfuAdjhE + 900*
      traceLambda12AdjfuTfuAdjLambda12 + 2700*traceYdAdjYuTYuAdjYd + 2700*
      traceYuAdjYdTYdAdjYu + 16200*traceYuAdjYuTYuAdjYu - 720*traceAdjYuTYu*Sqr
      (g1) + 720*MassB*traceYuAdjYu*Sqr(g1) - 1350*traceAdjfuTfu*Sqr(g1p) + 270
      *traceAdjYuTYu*Sqr(g1p) + 1350*MassBp*tracefuAdjfu*Sqr(g1p) - 270*MassBp*
      traceYuAdjYu*Sqr(g1p) + 483*MassB*Sqr(g1)*Sqr(g1p) + 483*MassBp*Sqr(g1)*
      Sqr(g1p) + 900*MassB*Sqr(g1)*Sqr(g2) + 900*MassWB*Sqr(g1)*Sqr(g2) + 675*
      MassBp*Sqr(g1p)*Sqr(g2) + 675*MassWB*Sqr(g1p)*Sqr(g2) - 14400*
      traceAdjYuTYu*Sqr(g3) + 14400*MassG*traceYuAdjYu*Sqr(g3) + 2720*MassB*Sqr
      (g1)*Sqr(g3) + 2720*MassG*Sqr(g1)*Sqr(g3) + 480*MassBp*Sqr(g1p)*Sqr(g3) +
      480*MassG*Sqr(g1p)*Sqr(g3) + 7200*MassG*Sqr(g2)*Sqr(g3) + 7200*MassWB*
      Sqr(g2)*Sqr(g3) + 27*Power(g1p,4)*MassBp*Sqr(QS) + 5400*Lambdax*Sqr(Conj(
      Lambdax))*TLambdax + 45*Conj(Lambdax)*((20*tracefdAdjfd + 60*
      traceKappaAdjKappa + 40*traceLambda12AdjLambda12 + 60*traceYdAdjYd + 20*
      traceYeAdjYe + 20*AbsSqr(Sigmax) - 5*Sqr(g1p) - Sqr(g1p)*Sqr(QS))*
      TLambdax + Lambdax*(20*traceAdjfdTfd + 60*traceAdjKappaTKappa + 40*
      traceAdjLambda12TLambda12 + 60*traceAdjYdTYd + 20*traceAdjYeTYe + 5*
      MassBp*Sqr(g1p) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) -
      0.4*(5*traceAdjfdTfd + 15*traceAdjYdTYd + 5*traceAdjYeTYe + 2*MassB*Sqr(
      g1) + 3*MassBp*Sqr(g1p) + 5*Conj(Lambdax)*TLambdax)*(Yu*Yd.adjoint()*Yd)
      - 2*tracefdAdjfd*(Yu*Yd.adjoint()*TYd) - 6*traceYdAdjYd*(Yu*Yd.adjoint()*
      TYd) - 2*traceYeAdjYe*(Yu*Yd.adjoint()*TYd) - 2*AbsSqr(Lambdax)*(Yu*
      Yd.adjoint()*TYd) + 0.8*Sqr(g1)*(Yu*Yd.adjoint()*TYd) + 1.2*Sqr(g1p)*(Yu*
      Yd.adjoint()*TYd) - 6*traceAdjfuTfu*(Yu*Yu.adjoint()*Yu) - 18*
      traceAdjYuTYu*(Yu*Yu.adjoint()*Yu) - 0.8*MassB*Sqr(g1)*(Yu*Yu.adjoint()*
      Yu) - 1.2*MassBp*Sqr(g1p)*(Yu*Yu.adjoint()*Yu) - 12*MassWB*Sqr(g2)*(Yu*
      Yu.adjoint()*Yu) - 6*Conj(Lambdax)*TLambdax*(Yu*Yu.adjoint()*Yu) - 4*
      tracefuAdjfu*(Yu*Yu.adjoint()*TYu) - 12*traceYuAdjYu*(Yu*Yu.adjoint()*TYu
      ) - 4*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*TYu) + 1.2*Sqr(g1)*(Yu*Yu.adjoint(
      )*TYu) + 0.8*Sqr(g1p)*(Yu*Yu.adjoint()*TYu) + 6*Sqr(g2)*(Yu*Yu.adjoint()*
      TYu) - 6*traceAdjgDTgD*(Yu*gD.conjugate()*gD.transpose()) - 2*
      traceAdjhEThE*(Yu*gD.conjugate()*gD.transpose()) - 0.8*MassB*Sqr(g1)*(Yu*
      gD.conjugate()*gD.transpose()) - 1.2*MassBp*Sqr(g1p)*(Yu*gD.conjugate()*
      gD.transpose()) - 2*Conj(SigmaL)*TSigmaL*(Yu*gD.conjugate()*gD.transpose(
      )) - 6*tracegDAdjgD*(Yu*gD.conjugate()*(TgD).transpose()) - 2*
      tracehEAdjhE*(Yu*gD.conjugate()*(TgD).transpose()) - 2*AbsSqr(SigmaL)*(Yu
      *gD.conjugate()*(TgD).transpose()) + 0.8*Sqr(g1)*(Yu*gD.conjugate()*(TgD)
      .transpose()) + 1.2*Sqr(g1p)*(Yu*gD.conjugate()*(TgD).transpose()) -
      tracefdAdjfd*(TYu*Yd.adjoint()*Yd) - 3*traceYdAdjYd*(TYu*Yd.adjoint()*Yd)
      - traceYeAdjYe*(TYu*Yd.adjoint()*Yd) - AbsSqr(Lambdax)*(TYu*Yd.adjoint()
      *Yd) + 0.4*Sqr(g1)*(TYu*Yd.adjoint()*Yd) + 0.6*Sqr(g1p)*(TYu*Yd.adjoint()
      *Yd) - 5*tracefuAdjfu*(TYu*Yu.adjoint()*Yu) - 15*traceYuAdjYu*(TYu*
      Yu.adjoint()*Yu) - 5*AbsSqr(Lambdax)*(TYu*Yu.adjoint()*Yu) + Sqr(g1p)*(
      TYu*Yu.adjoint()*Yu) + 12*Sqr(g2)*(TYu*Yu.adjoint()*Yu) - 3*tracegDAdjgD*
      (TYu*gD.conjugate()*gD.transpose()) - tracehEAdjhE*(TYu*gD.conjugate()*
      gD.transpose()) - AbsSqr(SigmaL)*(TYu*gD.conjugate()*gD.transpose()) +
      0.4*Sqr(g1)*(TYu*gD.conjugate()*gD.transpose()) + 0.6*Sqr(g1p)*(TYu*
      gD.conjugate()*gD.transpose()) - 4*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd)
      - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 2*(Yu*gD.conjugate()*gD.transpose()*Yu.adjoint()*TYu) - 4*(Yu*
      gD.conjugate()*gD.transpose()*gD.conjugate()*(TgD).transpose()) - 2*(Yu*
      gD.conjugate()*(Kappa).transpose()*Kappa.conjugate()*(TgD).transpose()) -
      4*(Yu*gD.conjugate()*(TgD).transpose()*Yu.adjoint()*Yu) - 4*(Yu*
      gD.conjugate()*(TgD).transpose()*gD.conjugate()*gD.transpose()) - 2*(Yu*
      gD.conjugate()*(TKappa).transpose()*Kappa.conjugate()*gD.transpose()) - 2
      *(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 4*(TYu*
      gD.conjugate()*gD.transpose()*Yu.adjoint()*Yu) - 2*(TYu*gD.conjugate()*
      gD.transpose()*gD.conjugate()*gD.transpose()) - TYu*gD.conjugate()*(Kappa
      ).transpose()*Kappa.conjugate()*gD.transpose())).real();


   return beta_TYu;
}

} // namespace flexiblesusy
