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

// File generated at Wed 3 Jun 2015 23:43:01

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(tracefdAdjfd*TYe + 3*traceYdAdjYd*TYe +
      traceYeAdjYe*TYe + AbsSqr(Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 0.7*Sqr(g1p)*
      TYe - 3*Sqr(g2)*TYe + Ye*(2*traceAdjfdTfd + 6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 3.6*MassB*Sqr(g1) + 1.4*MassBp*Sqr(g1p) + 6*MassWB*Sqr(g2
      ) + 2*Conj(Lambdax)*TLambdax) + 2*(hE*hE.adjoint()*TYe) + 4*(Ye*
      Ye.adjoint()*TYe) + 4*(ThE*hE.adjoint()*Ye) + 5*(TYe*Ye.adjoint()*Ye)))
      .real();


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(18.9*Power(g1,4)*TYe + 6.825*Power(g1p,4)*TYe +
      16.5*Power(g2,4)*TYe - 3*tracefdAdjfdfdAdjfd*TYe - 2*tracefdAdjfdfuAdjfu*
      TYe - 3*tracegDAdjgDTpYdconjYd*TYe - 2*tracehEAdjhEYeAdjYe*TYe -
      traceLambda12AdjLambda12Tpfdconjfd*TYe - 9*traceYdAdjYdYdAdjYd*TYe - 3*
      traceYdAdjYuYuAdjYd*TYe - 3*traceYeAdjYeYeAdjYe*TYe - tracefuAdjfu*AbsSqr
      (Lambdax)*TYe - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*TYe - 2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYe - 3*traceYuAdjYu*AbsSqr(
      Lambdax)*TYe - AbsSqr(Lambdax)*AbsSqr(Sigmax)*TYe - 0.4*traceYdAdjYd*Sqr(
      g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe + tracefdAdjfd*Sqr(g1p)*TYe - 0.6*
      traceYdAdjYd*Sqr(g1p)*TYe - 0.2*traceYeAdjYe*Sqr(g1p)*TYe - 0.25*AbsSqr(
      Lambdax)*Sqr(g1p)*TYe + 0.15*Sqr(g1)*Sqr(g1p)*TYe + 1.8*Sqr(g1)*Sqr(g2)*
      TYe + 1.95*Sqr(g1p)*Sqr(g2)*TYe + 16*traceYdAdjYd*Sqr(g3)*TYe + 0.035*
      Power(g1p,4)*Sqr(QS)*TYe + 0.05*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TYe - 3*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe - 0.02*Ye*(3780*Power(g1,4)*MassB +
      1365*Power(g1p,4)*MassBp + 3300*Power(g2,4)*MassWB + 100*
      traceAdjfdTfdconjLambda12TpLambda12 + 300*traceAdjgDTpYdconjYdTgD + 100*
      traceAdjLambda12TpfdconjfdTLambda12 + 300*traceAdjYdTYdconjgDTpgD + 600*
      tracefdAdjfdTfdAdjfd + 200*tracefdAdjfdTfuAdjfu + 200*
      tracefuAdjfuTfdAdjfd + 200*tracehEAdjhETYeAdjYe + 1800*
      traceYdAdjYdTYdAdjYd + 300*traceYdAdjYuTYuAdjYd + 200*
      traceYeAdjYeThEAdjhE + 600*traceYeAdjYeTYeAdjYe + 300*
      traceYuAdjYdTYdAdjYu + 40*traceAdjYdTYd*Sqr(g1) - 120*traceAdjYeTYe*Sqr(
      g1) - 40*MassB*traceYdAdjYd*Sqr(g1) + 120*MassB*traceYeAdjYe*Sqr(g1) -
      100*traceAdjfdTfd*Sqr(g1p) + 60*traceAdjYdTYd*Sqr(g1p) + 20*traceAdjYeTYe
      *Sqr(g1p) + 100*MassBp*tracefdAdjfd*Sqr(g1p) - 60*MassBp*traceYdAdjYd*Sqr
      (g1p) - 20*MassBp*traceYeAdjYe*Sqr(g1p) + 15*MassB*Sqr(g1)*Sqr(g1p) + 15*
      MassBp*Sqr(g1)*Sqr(g1p) + 180*MassB*Sqr(g1)*Sqr(g2) + 180*MassWB*Sqr(g1)*
      Sqr(g2) + 195*MassBp*Sqr(g1p)*Sqr(g2) + 195*MassWB*Sqr(g1p)*Sqr(g2) -
      1600*traceAdjYdTYd*Sqr(g3) + 1600*MassG*traceYdAdjYd*Sqr(g3) + 7*Power(
      g1p,4)*MassBp*Sqr(QS) + 600*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 5*Conj(
      Lambdax)*((20*tracefuAdjfu + 60*traceKappaAdjKappa + 40*
      traceLambda12AdjLambda12 + 60*traceYuAdjYu + 20*AbsSqr(Sigmax) + 5*Sqr(
      g1p) - Sqr(g1p)*Sqr(QS))*TLambdax + Lambdax*(20*traceAdjfuTfu + 60*
      traceAdjKappaTKappa + 40*traceAdjLambda12TLambda12 + 60*traceAdjYuTYu - 5
      *MassBp*Sqr(g1p) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) +
      0.8*(-15*traceAdjgDTgD - 5*traceAdjhEThE + 3*MassB*Sqr(g1) - 3*MassBp*
      Sqr(g1p) - 15*MassWB*Sqr(g2) - 5*Conj(SigmaL)*TSigmaL)*(hE*hE.adjoint()*
      Ye) - 6*tracegDAdjgD*(hE*hE.adjoint()*TYe) - 2*tracehEAdjhE*(hE*
      hE.adjoint()*TYe) - 2*AbsSqr(SigmaL)*(hE*hE.adjoint()*TYe) - 1.2*Sqr(g1)*
      (hE*hE.adjoint()*TYe) + 1.2*Sqr(g1p)*(hE*hE.adjoint()*TYe) + 6*Sqr(g2)*(
      hE*hE.adjoint()*TYe) - 6*traceAdjfdTfd*(Ye*Ye.adjoint()*Ye) - 18*
      traceAdjYdTYd*(Ye*Ye.adjoint()*Ye) - 6*traceAdjYeTYe*(Ye*Ye.adjoint()*Ye)
      - 3*MassBp*Sqr(g1p)*(Ye*Ye.adjoint()*Ye) - 12*MassWB*Sqr(g2)*(Ye*
      Ye.adjoint()*Ye) - 6*Conj(Lambdax)*TLambdax*(Ye*Ye.adjoint()*Ye) - 4*
      tracefdAdjfd*(Ye*Ye.adjoint()*TYe) - 12*traceYdAdjYd*(Ye*Ye.adjoint()*TYe
      ) - 4*traceYeAdjYe*(Ye*Ye.adjoint()*TYe) - 4*AbsSqr(Lambdax)*(Ye*
      Ye.adjoint()*TYe) + 1.2*Sqr(g1)*(Ye*Ye.adjoint()*TYe) + 1.8*Sqr(g1p)*(Ye*
      Ye.adjoint()*TYe) + 6*Sqr(g2)*(Ye*Ye.adjoint()*TYe) - 12*tracegDAdjgD*(
      ThE*hE.adjoint()*Ye) - 4*tracehEAdjhE*(ThE*hE.adjoint()*Ye) - 4*AbsSqr(
      SigmaL)*(ThE*hE.adjoint()*Ye) - 2.4*Sqr(g1)*(ThE*hE.adjoint()*Ye) + 2.4*
      Sqr(g1p)*(ThE*hE.adjoint()*Ye) + 12*Sqr(g2)*(ThE*hE.adjoint()*Ye) - 5*
      tracefdAdjfd*(TYe*Ye.adjoint()*Ye) - 15*traceYdAdjYd*(TYe*Ye.adjoint()*Ye
      ) - 5*traceYeAdjYe*(TYe*Ye.adjoint()*Ye) - 5*AbsSqr(Lambdax)*(TYe*
      Ye.adjoint()*Ye) - 1.2*Sqr(g1)*(TYe*Ye.adjoint()*Ye) + 2.7*Sqr(g1p)*(TYe*
      Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()*Ye) - 2*(hE*fu.adjoint()*
      fu*hE.adjoint()*TYe) - 4*(hE*fu.adjoint()*Tfu*hE.adjoint()*Ye) - 2*(hE*
      hE.adjoint()*hE*hE.adjoint()*TYe) - 4*(hE*hE.adjoint()*ThE*hE.adjoint()*
      Ye) - 2*(hE*(Lambda12).adjoint()*Lambda12*hE.adjoint()*TYe) - 4*(hE*(
      Lambda12).adjoint()*TLambda12*hE.adjoint()*Ye) - 4*(Ye*Ye.adjoint()*hE*
      hE.adjoint()*TYe) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 4*(Ye*
      Ye.adjoint()*ThE*hE.adjoint()*Ye) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*
      Ye) - 4*(ThE*fu.adjoint()*fu*hE.adjoint()*Ye) - 4*(ThE*hE.adjoint()*hE*
      hE.adjoint()*Ye) - 4*(ThE*(Lambda12).adjoint()*Lambda12*hE.adjoint()*Ye)
      - 2*(TYe*Ye.adjoint()*hE*hE.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

} // namespace flexiblesusy
