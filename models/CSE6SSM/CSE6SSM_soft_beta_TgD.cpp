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

// File generated at Wed 3 Jun 2015 23:43:05

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TgD.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TgD_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,3> beta_TgD;

   beta_TgD = (oneOver16PiSqr*(3*tracegDAdjgD*TgD + tracehEAdjhE*TgD +
      AbsSqr(SigmaL)*TgD - 0.4666666666666667*Sqr(g1)*TgD - 0.7*Sqr(g1p)*TgD -
      3*Sqr(g2)*TgD - 5.333333333333333*Sqr(g3)*TgD + gD*(6*traceAdjgDTgD + 2*
      traceAdjhEThE + 0.9333333333333333*MassB*Sqr(g1) + 1.4*MassBp*Sqr(g1p) +
      6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*Conj(SigmaL)*
      TSigmaL) + 5*(gD*gD.adjoint()*TgD) + 2*(gD*(Kappa).adjoint()*TKappa) + 4*
      (TgD*gD.adjoint()*gD) + TgD*(Kappa).adjoint()*Kappa + Yd.transpose()*
      Yd.conjugate()*TgD + Yu.transpose()*Yu.conjugate()*TgD + 2*((TYd)
      .transpose()*Yd.conjugate()*gD) + 2*((TYu).transpose()*Yu.conjugate()*gD)
      )).real();


   return beta_TgD;
}

/**
 * Calculates the two-loop beta function of TgD.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TgD_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracefuAdjhEThEAdjfu = TRACE_STRUCT.tracefuAdjhEThEAdjfu;
   const double tracegDAdjgDTgDAdjgD = TRACE_STRUCT.tracegDAdjgDTgDAdjgD;
   const double tracegDAdjKappaTKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD;
   const double tracehEAdjfuTfuAdjhE = TRACE_STRUCT.tracehEAdjfuTfuAdjhE;
   const double tracehEAdjhEThEAdjhE = TRACE_STRUCT.tracehEAdjhEThEAdjhE;
   const double tracehEAdjhETYeAdjYe = TRACE_STRUCT.tracehEAdjhETYeAdjYe;
   const double tracehEAdjLambda12TLambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE;
   const double traceYeAdjYeThEAdjhE = TRACE_STRUCT.traceYeAdjYeThEAdjhE;
   const double traceKappaAdjgDTgDAdjKappa =
      TRACE_STRUCT.traceKappaAdjgDTgDAdjKappa;
   const double traceLambda12AdjhEThEAdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjhEThEAdjLambda12;
   const double traceAdjgDTpYdconjYdTgD =
      TRACE_STRUCT.traceAdjgDTpYdconjYdTgD;
   const double traceAdjgDTpYuconjYuTgD =
      TRACE_STRUCT.traceAdjgDTpYuconjYuTgD;
   const double traceAdjYdTYdconjgDTpgD =
      TRACE_STRUCT.traceAdjYdTYdconjgDTpgD;
   const double traceAdjYuTYuconjgDTpgD =
      TRACE_STRUCT.traceAdjYuTYuconjgDTpgD;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracefuAdjhEhEAdjfu = TRACE_STRUCT.tracefuAdjhEhEAdjfu;
   const double tracegDAdjgDgDAdjgD = TRACE_STRUCT.tracegDAdjgDgDAdjgD;
   const double tracegDAdjgDTpYdconjYd =
      TRACE_STRUCT.tracegDAdjgDTpYdconjYd;
   const double tracegDAdjgDTpYuconjYu =
      TRACE_STRUCT.tracegDAdjgDTpYuconjYu;
   const double tracegDAdjKappaKappaAdjgD =
      TRACE_STRUCT.tracegDAdjKappaKappaAdjgD;
   const double tracehEAdjhEhEAdjhE = TRACE_STRUCT.tracehEAdjhEhEAdjhE;
   const double tracehEAdjhEYeAdjYe = TRACE_STRUCT.tracehEAdjhEYeAdjYe;
   const double tracehEAdjLambda12Lambda12AdjhE =
      TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE;


   Eigen::Matrix<double,3,3> beta_TgD;

   beta_TgD = (twoLoop*(4.588888888888889*Power(g1,4)*TgD + 6.825*Power(
      g1p,4)*TgD + 16.5*Power(g2,4)*TgD + 14.222222222222221*Power(g3,4)*TgD -
      tracefuAdjhEhEAdjfu*TgD - 9*tracegDAdjgDgDAdjgD*TgD - 3*
      tracegDAdjgDTpYdconjYd*TgD - 3*tracegDAdjgDTpYuconjYu*TgD - 3*
      tracegDAdjKappaKappaAdjgD*TgD - 3*tracehEAdjhEhEAdjhE*TgD - 2*
      tracehEAdjhEYeAdjYe*TgD - tracehEAdjLambda12Lambda12AdjhE*TgD - 2*AbsSqr(
      KappaPr)*AbsSqr(SigmaL)*TgD - AbsSqr(Sigmax)*AbsSqr(SigmaL)*TgD - 0.4*
      tracegDAdjgD*Sqr(g1)*TgD + 1.2*tracehEAdjhE*Sqr(g1)*TgD + 0.9*
      tracegDAdjgD*Sqr(g1p)*TgD + 0.3*tracehEAdjhE*Sqr(g1p)*TgD +
      0.6833333333333333*Sqr(g1)*Sqr(g1p)*TgD + Sqr(g1)*Sqr(g2)*TgD + 0.75*Sqr(
      g1p)*Sqr(g2)*TgD + 16*tracegDAdjgD*Sqr(g3)*TgD + 0.8888888888888888*Sqr(
      g1)*Sqr(g3)*TgD + 2.6666666666666665*Sqr(g1p)*Sqr(g3)*TgD + 8*Sqr(g2)*Sqr
      (g3)*TgD + 0.035*Power(g1p,4)*Sqr(QS)*TgD - 3*Sqr(Conj(SigmaL))*Sqr(
      SigmaL)*TgD + gD*(-18.355555555555554*Power(g1,4)*MassB - 27.3*Power(g1p,
      4)*MassBp - 56.888888888888886*Power(g3,4)*MassG - 66*Power(g2,4)*MassWB
      - 6*traceAdjgDTpYdconjYdTgD - 6*traceAdjgDTpYuconjYuTgD - 6*
      traceAdjYdTYdconjgDTpgD - 6*traceAdjYuTYuconjgDTpgD - 2*
      tracefuAdjhEThEAdjfu - 36*tracegDAdjgDTgDAdjgD - 6*
      tracegDAdjKappaTKappaAdjgD - 2*tracehEAdjfuTfuAdjhE - 12*
      tracehEAdjhEThEAdjhE - 4*tracehEAdjhETYeAdjYe - 2*
      tracehEAdjLambda12TLambda12AdjhE - 6*traceKappaAdjgDTgDAdjKappa - 2*
      traceLambda12AdjhEThEAdjLambda12 - 4*traceYeAdjYeThEAdjhE - 0.8*
      traceAdjgDTgD*Sqr(g1) + 2.4*traceAdjhEThE*Sqr(g1) + 0.8*MassB*
      tracegDAdjgD*Sqr(g1) - 2.4*MassB*tracehEAdjhE*Sqr(g1) + 1.8*traceAdjgDTgD
      *Sqr(g1p) + 0.6*traceAdjhEThE*Sqr(g1p) - 1.8*MassBp*tracegDAdjgD*Sqr(g1p)
      - 0.6*MassBp*tracehEAdjhE*Sqr(g1p) - 1.3666666666666667*MassB*Sqr(g1)*
      Sqr(g1p) - 1.3666666666666667*MassBp*Sqr(g1)*Sqr(g1p) - 2*MassB*Sqr(g1)*
      Sqr(g2) - 2*MassWB*Sqr(g1)*Sqr(g2) - 1.5*MassBp*Sqr(g1p)*Sqr(g2) - 1.5*
      MassWB*Sqr(g1p)*Sqr(g2) + 32*traceAdjgDTgD*Sqr(g3) - 32*MassG*
      tracegDAdjgD*Sqr(g3) - 1.7777777777777777*MassB*Sqr(g1)*Sqr(g3) -
      1.7777777777777777*MassG*Sqr(g1)*Sqr(g3) - 5.333333333333333*MassBp*Sqr(
      g1p)*Sqr(g3) - 5.333333333333333*MassG*Sqr(g1p)*Sqr(g3) - 16*MassG*Sqr(g2
      )*Sqr(g3) - 16*MassWB*Sqr(g2)*Sqr(g3) - 0.14*Power(g1p,4)*MassBp*Sqr(QS)
      - 12*SigmaL*Sqr(Conj(SigmaL))*TSigmaL - 4*Conj(KappaPr)*Conj(SigmaL)*(
      SigmaL*TKappaPr + KappaPr*TSigmaL) - 2*Conj(Sigmax)*Conj(SigmaL)*(SigmaL*
      TSigmax + Sigmax*TSigmaL)) - 0.4*(45*traceAdjgDTgD + 15*traceAdjhEThE + 4
      *MassB*Sqr(g1) + MassBp*Sqr(g1p) + 30*MassWB*Sqr(g2) + 15*Conj(SigmaL)*
      TSigmaL)*(gD*gD.adjoint()*gD) - 15*tracegDAdjgD*(gD*gD.adjoint()*TgD) - 5
      *tracehEAdjhE*(gD*gD.adjoint()*TgD) - 5*AbsSqr(SigmaL)*(gD*gD.adjoint()*
      TgD) + 1.2*Sqr(g1)*(gD*gD.adjoint()*TgD) - 0.2*Sqr(g1p)*(gD*gD.adjoint()*
      TgD) + 12*Sqr(g2)*(gD*gD.adjoint()*TgD) - 6*traceAdjKappaTKappa*(gD*(
      Kappa).adjoint()*Kappa) - 4*traceAdjLambda12TLambda12*(gD*(Kappa).adjoint
      ()*Kappa) + 0.5*MassBp*Sqr(g1p)*(gD*(Kappa).adjoint()*Kappa) - 0.1*MassBp
      *Sqr(g1p)*Sqr(QS)*(gD*(Kappa).adjoint()*Kappa) - 4*Conj(Lambdax)*TLambdax
      *(gD*(Kappa).adjoint()*Kappa) - 2*Conj(Sigmax)*TSigmax*(gD*(Kappa)
      .adjoint()*Kappa) - 6*traceKappaAdjKappa*(gD*(Kappa).adjoint()*TKappa) -
      4*traceLambda12AdjLambda12*(gD*(Kappa).adjoint()*TKappa) - 4*AbsSqr(
      Lambdax)*(gD*(Kappa).adjoint()*TKappa) - 2*AbsSqr(Sigmax)*(gD*(Kappa)
      .adjoint()*TKappa) - 0.5*Sqr(g1p)*(gD*(Kappa).adjoint()*TKappa) + 0.1*Sqr
      (g1p)*Sqr(QS)*(gD*(Kappa).adjoint()*TKappa) - 12*tracegDAdjgD*(TgD*
      gD.adjoint()*gD) - 4*tracehEAdjhE*(TgD*gD.adjoint()*gD) - 4*AbsSqr(SigmaL
      )*(TgD*gD.adjoint()*gD) + 1.2*Sqr(g1)*(TgD*gD.adjoint()*gD) + 0.8*Sqr(g1p
      )*(TgD*gD.adjoint()*gD) + 6*Sqr(g2)*(TgD*gD.adjoint()*gD) - 3*
      traceKappaAdjKappa*(TgD*(Kappa).adjoint()*Kappa) - 2*
      traceLambda12AdjLambda12*(TgD*(Kappa).adjoint()*Kappa) - 2*AbsSqr(Lambdax
      )*(TgD*(Kappa).adjoint()*Kappa) - AbsSqr(Sigmax)*(TgD*(Kappa).adjoint()*
      Kappa) - 0.25*Sqr(g1p)*(TgD*(Kappa).adjoint()*Kappa) + 0.05*Sqr(g1p)*Sqr(
      QS)*(TgD*(Kappa).adjoint()*Kappa) - 2*traceAdjfdTfd*(Yd.transpose()*
      Yd.conjugate()*gD) - 6*traceAdjYdTYd*(Yd.transpose()*Yd.conjugate()*gD) -
      2*traceAdjYeTYe*(Yd.transpose()*Yd.conjugate()*gD) - 0.8*MassB*Sqr(g1)*(
      Yd.transpose()*Yd.conjugate()*gD) - 1.2*MassBp*Sqr(g1p)*(Yd.transpose()*
      Yd.conjugate()*gD) - 2*Conj(Lambdax)*TLambdax*(Yd.transpose()*
      Yd.conjugate()*gD) - tracefdAdjfd*(Yd.transpose()*Yd.conjugate()*TgD) - 3
      *traceYdAdjYd*(Yd.transpose()*Yd.conjugate()*TgD) - traceYeAdjYe*(
      Yd.transpose()*Yd.conjugate()*TgD) - AbsSqr(Lambdax)*(Yd.transpose()*
      Yd.conjugate()*TgD) + 0.4*Sqr(g1)*(Yd.transpose()*Yd.conjugate()*TgD) +
      0.6*Sqr(g1p)*(Yd.transpose()*Yd.conjugate()*TgD) - 2*traceAdjfuTfu*(
      Yu.transpose()*Yu.conjugate()*gD) - 6*traceAdjYuTYu*(Yu.transpose()*
      Yu.conjugate()*gD) - 1.6*MassB*Sqr(g1)*(Yu.transpose()*Yu.conjugate()*gD)
      - 0.4*MassBp*Sqr(g1p)*(Yu.transpose()*Yu.conjugate()*gD) - 2*Conj(
      Lambdax)*TLambdax*(Yu.transpose()*Yu.conjugate()*gD) - tracefuAdjfu*(
      Yu.transpose()*Yu.conjugate()*TgD) - 3*traceYuAdjYu*(Yu.transpose()*
      Yu.conjugate()*TgD) - AbsSqr(Lambdax)*(Yu.transpose()*Yu.conjugate()*TgD)
      + 0.8*Sqr(g1)*(Yu.transpose()*Yu.conjugate()*TgD) + 0.2*Sqr(g1p)*(
      Yu.transpose()*Yu.conjugate()*TgD) - 2*tracefdAdjfd*((TYd).transpose()*
      Yd.conjugate()*gD) - 6*traceYdAdjYd*((TYd).transpose()*Yd.conjugate()*gD)
      - 2*traceYeAdjYe*((TYd).transpose()*Yd.conjugate()*gD) - 2*AbsSqr(
      Lambdax)*((TYd).transpose()*Yd.conjugate()*gD) + 0.8*Sqr(g1)*((TYd)
      .transpose()*Yd.conjugate()*gD) + 1.2*Sqr(g1p)*((TYd).transpose()*
      Yd.conjugate()*gD) - 2*tracefuAdjfu*((TYu).transpose()*Yu.conjugate()*gD)
      - 6*traceYuAdjYu*((TYu).transpose()*Yu.conjugate()*gD) - 2*AbsSqr(
      Lambdax)*((TYu).transpose()*Yu.conjugate()*gD) + 1.6*Sqr(g1)*((TYu)
      .transpose()*Yu.conjugate()*gD) + 0.4*Sqr(g1p)*((TYu).transpose()*
      Yu.conjugate()*gD) - 6*(gD*gD.adjoint()*gD*gD.adjoint()*TgD) - 8*(gD*
      gD.adjoint()*TgD*gD.adjoint()*gD) - 4*(gD*gD.adjoint()*Yd.transpose()*
      Yd.conjugate()*TgD) - 4*(gD*gD.adjoint()*Yu.transpose()*Yu.conjugate()*
      TgD) - 4*(gD*gD.adjoint()*(TYd).transpose()*Yd.conjugate()*gD) - 4*(gD*
      gD.adjoint()*(TYu).transpose()*Yu.conjugate()*gD) - gD*(Kappa).adjoint()*
      Kappa*gD.adjoint()*TgD - 2*(gD*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*
      TKappa) - 2*(gD*(Kappa).adjoint()*TKappa*gD.adjoint()*gD) - 2*(gD*(Kappa)
      .adjoint()*TKappa*(Kappa).adjoint()*Kappa) - 6*(TgD*gD.adjoint()*gD*
      gD.adjoint()*gD) - 2*(TgD*gD.adjoint()*Yd.transpose()*Yd.conjugate()*gD)
      - 2*(TgD*gD.adjoint()*Yu.transpose()*Yu.conjugate()*gD) - 2*(TgD*(Kappa)
      .adjoint()*Kappa*gD.adjoint()*gD) - TgD*(Kappa).adjoint()*Kappa*(Kappa)
      .adjoint()*Kappa - 2*(Yd.transpose()*Yd.conjugate()*Yd.transpose()*
      Yd.conjugate()*TgD) - 4*(Yd.transpose()*Yd.conjugate()*(TYd).transpose()*
      Yd.conjugate()*gD) - 2*(Yu.transpose()*Yu.conjugate()*Yu.transpose()*
      Yu.conjugate()*TgD) - 4*(Yu.transpose()*Yu.conjugate()*(TYu).transpose()*
      Yu.conjugate()*gD) - 4*((TYd).transpose()*Yd.conjugate()*Yd.transpose()*
      Yd.conjugate()*gD) - 4*((TYu).transpose()*Yu.conjugate()*Yu.transpose()*
      Yu.conjugate()*gD))).real();


   return beta_TgD;
}

} // namespace flexiblesusy
