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

// File generated at Wed 3 Jun 2015 23:42:59

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ThE.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_ThE_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;


   Eigen::Matrix<double,3,2> beta_ThE;

   beta_ThE = (oneOver16PiSqr*(3*tracegDAdjgD*ThE + tracehEAdjhE*ThE +
      AbsSqr(SigmaL)*ThE - 1.8*Sqr(g1)*ThE - 0.7*Sqr(g1p)*ThE - 3*Sqr(g2)*ThE +
      hE*(6*traceAdjgDTgD + 2*traceAdjhEThE + 3.6*MassB*Sqr(g1) + 1.4*MassBp*
      Sqr(g1p) + 6*MassWB*Sqr(g2) + 2*Conj(SigmaL)*TSigmaL) + 2*(hE*fu.adjoint(
      )*Tfu) + 4*(hE*hE.adjoint()*ThE) + 2*(hE*(Lambda12).adjoint()*TLambda12)
      + 2*(Ye*Ye.adjoint()*ThE) + ThE*fu.adjoint()*fu + 5*(ThE*hE.adjoint()*hE)
      + ThE*(Lambda12).adjoint()*Lambda12 + 4*(TYe*Ye.adjoint()*hE))).real();


   return beta_ThE;
}

/**
 * Calculates the two-loop beta function of ThE.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_ThE_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
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


   Eigen::Matrix<double,3,2> beta_ThE;

   beta_ThE = (twoLoop*(18.9*Power(g1,4)*ThE + 6.825*Power(g1p,4)*ThE +
      16.5*Power(g2,4)*ThE - tracefuAdjhEhEAdjfu*ThE - 9*tracegDAdjgDgDAdjgD*
      ThE - 3*tracegDAdjgDTpYdconjYd*ThE - 3*tracegDAdjgDTpYuconjYu*ThE - 3*
      tracegDAdjKappaKappaAdjgD*ThE - 3*tracehEAdjhEhEAdjhE*ThE - 2*
      tracehEAdjhEYeAdjYe*ThE - tracehEAdjLambda12Lambda12AdjhE*ThE - 2*AbsSqr(
      KappaPr)*AbsSqr(SigmaL)*ThE - AbsSqr(Sigmax)*AbsSqr(SigmaL)*ThE - 0.4*
      tracegDAdjgD*Sqr(g1)*ThE + 1.2*tracehEAdjhE*Sqr(g1)*ThE + 0.9*
      tracegDAdjgD*Sqr(g1p)*ThE + 0.3*tracehEAdjhE*Sqr(g1p)*ThE + 0.15*Sqr(g1)*
      Sqr(g1p)*ThE + 1.8*Sqr(g1)*Sqr(g2)*ThE + 1.95*Sqr(g1p)*Sqr(g2)*ThE + 16*
      tracegDAdjgD*Sqr(g3)*ThE + 0.035*Power(g1p,4)*Sqr(QS)*ThE - 3*Sqr(Conj(
      SigmaL))*Sqr(SigmaL)*ThE - 0.02*hE*(3780*Power(g1,4)*MassB + 1365*Power(
      g1p,4)*MassBp + 3300*Power(g2,4)*MassWB + 300*traceAdjgDTpYdconjYdTgD +
      300*traceAdjgDTpYuconjYuTgD + 300*traceAdjYdTYdconjgDTpgD + 300*
      traceAdjYuTYuconjgDTpgD + 100*tracefuAdjhEThEAdjfu + 1800*
      tracegDAdjgDTgDAdjgD + 300*tracegDAdjKappaTKappaAdjgD + 100*
      tracehEAdjfuTfuAdjhE + 600*tracehEAdjhEThEAdjhE + 200*
      tracehEAdjhETYeAdjYe + 100*tracehEAdjLambda12TLambda12AdjhE + 300*
      traceKappaAdjgDTgDAdjKappa + 100*traceLambda12AdjhEThEAdjLambda12 + 200*
      traceYeAdjYeThEAdjhE + 40*traceAdjgDTgD*Sqr(g1) - 120*traceAdjhEThE*Sqr(
      g1) - 40*MassB*tracegDAdjgD*Sqr(g1) + 120*MassB*tracehEAdjhE*Sqr(g1) - 90
      *traceAdjgDTgD*Sqr(g1p) - 30*traceAdjhEThE*Sqr(g1p) + 90*MassBp*
      tracegDAdjgD*Sqr(g1p) + 30*MassBp*tracehEAdjhE*Sqr(g1p) + 15*MassB*Sqr(g1
      )*Sqr(g1p) + 15*MassBp*Sqr(g1)*Sqr(g1p) + 180*MassB*Sqr(g1)*Sqr(g2) + 180
      *MassWB*Sqr(g1)*Sqr(g2) + 195*MassBp*Sqr(g1p)*Sqr(g2) + 195*MassWB*Sqr(
      g1p)*Sqr(g2) - 1600*traceAdjgDTgD*Sqr(g3) + 1600*MassG*tracegDAdjgD*Sqr(
      g3) + 7*Power(g1p,4)*MassBp*Sqr(QS) + 600*SigmaL*Sqr(Conj(SigmaL))*
      TSigmaL + 200*Conj(KappaPr)*Conj(SigmaL)*(SigmaL*TKappaPr + KappaPr*
      TSigmaL) + 100*Conj(Sigmax)*Conj(SigmaL)*(SigmaL*TSigmax + Sigmax*TSigmaL
      )) - 2*(traceAdjfuTfu + 3*traceAdjYuTYu + MassBp*Sqr(g1p) + Conj(Lambdax)
      *TLambdax)*(hE*fu.adjoint()*fu) - 2*tracefuAdjfu*(hE*fu.adjoint()*Tfu) -
      6*traceYuAdjYu*(hE*fu.adjoint()*Tfu) - 2*AbsSqr(Lambdax)*(hE*fu.adjoint()
      *Tfu) + 2*Sqr(g1p)*(hE*fu.adjoint()*Tfu) - 18*traceAdjgDTgD*(hE*
      hE.adjoint()*hE) - 6*traceAdjhEThE*(hE*hE.adjoint()*hE) - 2*MassBp*Sqr(
      g1p)*(hE*hE.adjoint()*hE) - 12*MassWB*Sqr(g2)*(hE*hE.adjoint()*hE) - 6*
      Conj(SigmaL)*TSigmaL*(hE*hE.adjoint()*hE) - 12*tracegDAdjgD*(hE*
      hE.adjoint()*ThE) - 4*tracehEAdjhE*(hE*hE.adjoint()*ThE) - 4*AbsSqr(
      SigmaL)*(hE*hE.adjoint()*ThE) + 1.2*Sqr(g1)*(hE*hE.adjoint()*ThE) + 0.8*
      Sqr(g1p)*(hE*hE.adjoint()*ThE) + 6*Sqr(g2)*(hE*hE.adjoint()*ThE) - 6*
      traceAdjKappaTKappa*(hE*(Lambda12).adjoint()*Lambda12) - 4*
      traceAdjLambda12TLambda12*(hE*(Lambda12).adjoint()*Lambda12) + 0.5*MassBp
      *Sqr(g1p)*(hE*(Lambda12).adjoint()*Lambda12) - 0.1*MassBp*Sqr(g1p)*Sqr(QS
      )*(hE*(Lambda12).adjoint()*Lambda12) - 4*Conj(Lambdax)*TLambdax*(hE*(
      Lambda12).adjoint()*Lambda12) - 2*Conj(Sigmax)*TSigmax*(hE*(Lambda12)
      .adjoint()*Lambda12) - 6*traceKappaAdjKappa*(hE*(Lambda12).adjoint()*
      TLambda12) - 4*traceLambda12AdjLambda12*(hE*(Lambda12).adjoint()*
      TLambda12) - 4*AbsSqr(Lambdax)*(hE*(Lambda12).adjoint()*TLambda12) - 2*
      AbsSqr(Sigmax)*(hE*(Lambda12).adjoint()*TLambda12) - 0.5*Sqr(g1p)*(hE*(
      Lambda12).adjoint()*TLambda12) + 0.1*Sqr(g1p)*Sqr(QS)*(hE*(Lambda12)
      .adjoint()*TLambda12) - 4*traceAdjfdTfd*(Ye*Ye.adjoint()*hE) - 12*
      traceAdjYdTYd*(Ye*Ye.adjoint()*hE) - 4*traceAdjYeTYe*(Ye*Ye.adjoint()*hE)
      + 2.4*MassB*Sqr(g1)*(Ye*Ye.adjoint()*hE) - 2.4*MassBp*Sqr(g1p)*(Ye*
      Ye.adjoint()*hE) - 12*MassWB*Sqr(g2)*(Ye*Ye.adjoint()*hE) - 4*Conj(
      Lambdax)*TLambdax*(Ye*Ye.adjoint()*hE) - 2*tracefdAdjfd*(Ye*Ye.adjoint()*
      ThE) - 6*traceYdAdjYd*(Ye*Ye.adjoint()*ThE) - 2*traceYeAdjYe*(Ye*
      Ye.adjoint()*ThE) - 2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()*ThE) - 1.2*Sqr(g1)
      *(Ye*Ye.adjoint()*ThE) + 1.2*Sqr(g1p)*(Ye*Ye.adjoint()*ThE) + 6*Sqr(g2)*(
      Ye*Ye.adjoint()*ThE) - tracefuAdjfu*(ThE*fu.adjoint()*fu) - 3*
      traceYuAdjYu*(ThE*fu.adjoint()*fu) - AbsSqr(Lambdax)*(ThE*fu.adjoint()*fu
      ) + Sqr(g1p)*(ThE*fu.adjoint()*fu) - 15*tracegDAdjgD*(ThE*hE.adjoint()*hE
      ) - 5*tracehEAdjhE*(ThE*hE.adjoint()*hE) - 5*AbsSqr(SigmaL)*(ThE*
      hE.adjoint()*hE) - 1.2*Sqr(g1)*(ThE*hE.adjoint()*hE) + 2.2*Sqr(g1p)*(ThE*
      hE.adjoint()*hE) + 12*Sqr(g2)*(ThE*hE.adjoint()*hE) - 3*
      traceKappaAdjKappa*(ThE*(Lambda12).adjoint()*Lambda12) - 2*
      traceLambda12AdjLambda12*(ThE*(Lambda12).adjoint()*Lambda12) - 2*AbsSqr(
      Lambdax)*(ThE*(Lambda12).adjoint()*Lambda12) - AbsSqr(Sigmax)*(ThE*(
      Lambda12).adjoint()*Lambda12) - 0.25*Sqr(g1p)*(ThE*(Lambda12).adjoint()*
      Lambda12) + 0.05*Sqr(g1p)*Sqr(QS)*(ThE*(Lambda12).adjoint()*Lambda12) - 4
      *tracefdAdjfd*(TYe*Ye.adjoint()*hE) - 12*traceYdAdjYd*(TYe*Ye.adjoint()*
      hE) - 4*traceYeAdjYe*(TYe*Ye.adjoint()*hE) - 4*AbsSqr(Lambdax)*(TYe*
      Ye.adjoint()*hE) - 2.4*Sqr(g1)*(TYe*Ye.adjoint()*hE) + 2.4*Sqr(g1p)*(TYe*
      Ye.adjoint()*hE) + 12*Sqr(g2)*(TYe*Ye.adjoint()*hE) - 4*(hE*fu.adjoint()*
      fd*fd.adjoint()*Tfu) - 4*(hE*fu.adjoint()*fu*fu.adjoint()*Tfu) - 2*(hE*
      fu.adjoint()*fu*hE.adjoint()*ThE) - 4*(hE*fu.adjoint()*Tfd*fd.adjoint()*
      fu) - 4*(hE*fu.adjoint()*Tfu*fu.adjoint()*fu) - 4*(hE*fu.adjoint()*Tfu*
      hE.adjoint()*hE) - 6*(hE*hE.adjoint()*hE*hE.adjoint()*ThE) - 4*(hE*
      hE.adjoint()*Ye*Ye.adjoint()*ThE) - 8*(hE*hE.adjoint()*ThE*hE.adjoint()*
      hE) - 4*(hE*hE.adjoint()*TYe*Ye.adjoint()*hE) - 2*(hE*(Lambda12).adjoint(
      )*Lambda12*hE.adjoint()*ThE) - 2*(hE*(Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*TLambda12) - 4*(hE*(Lambda12).adjoint()*TLambda12*
      hE.adjoint()*hE) - 2*(hE*(Lambda12).adjoint()*TLambda12*(Lambda12)
      .adjoint()*Lambda12) - 2*(hE*(Lambda12).adjoint()*fd.transpose()*
      fd.conjugate()*TLambda12) - 2*(hE*(Lambda12).adjoint()*(Tfd).transpose()*
      fd.conjugate()*Lambda12) - 2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*ThE) - 4*(
      Ye*Ye.adjoint()*TYe*Ye.adjoint()*hE) - 2*(ThE*fu.adjoint()*fd*fd.adjoint(
      )*fu) - 2*(ThE*fu.adjoint()*fu*fu.adjoint()*fu) - 4*(ThE*fu.adjoint()*fu*
      hE.adjoint()*hE) - 6*(ThE*hE.adjoint()*hE*hE.adjoint()*hE) - 2*(ThE*
      hE.adjoint()*Ye*Ye.adjoint()*hE) - 4*(ThE*(Lambda12).adjoint()*Lambda12*
      hE.adjoint()*hE) - ThE*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()
      *Lambda12 - ThE*(Lambda12).adjoint()*fd.transpose()*fd.conjugate()*
      Lambda12 - 4*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*hE))).real();


   return beta_ThE;
}

} // namespace flexiblesusy
