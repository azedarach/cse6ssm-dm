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

// File generated at Wed 3 Jun 2015 23:43:13

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of Tfd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_Tfd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,2> beta_Tfd;

   beta_Tfd = (oneOver16PiSqr*(tracefdAdjfd*Tfd + 3*traceYdAdjYd*Tfd +
      traceYeAdjYe*Tfd + AbsSqr(Lambdax)*Tfd - 0.6*Sqr(g1)*Tfd - 1.9*Sqr(g1p)*
      Tfd - 3*Sqr(g2)*Tfd + fd*(2*traceAdjfdTfd + 6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 1.2*MassB*Sqr(g1) + 3.8*MassBp*Sqr(g1p) + 6*MassWB*Sqr(g2
      ) + 2*Conj(Lambdax)*TLambdax) + 4*(fd*fd.adjoint()*Tfd) + 2*(fd*
      Lambda12.conjugate()*(TLambda12).transpose()) + 2*(fu*fu.adjoint()*Tfd) +
      5*(Tfd*fd.adjoint()*fd) + Tfd*Lambda12.conjugate()*(Lambda12).transpose(
      ) + 4*(Tfu*fu.adjoint()*fd))).real();


   return beta_Tfd;
}

/**
 * Calculates the two-loop beta function of Tfd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,2> CSE6SSM_soft_parameters::calc_beta_Tfd_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,2> beta_Tfd;

   beta_Tfd = (twoLoop*(5.94*Power(g1,4)*Tfd + 19.665*Power(g1p,4)*Tfd +
      16.5*Power(g2,4)*Tfd - 3*tracefdAdjfdfdAdjfd*Tfd - 2*tracefdAdjfdfuAdjfu*
      Tfd - 3*tracegDAdjgDTpYdconjYd*Tfd - 2*tracehEAdjhEYeAdjYe*Tfd -
      traceLambda12AdjLambda12Tpfdconjfd*Tfd - 9*traceYdAdjYdYdAdjYd*Tfd - 3*
      traceYdAdjYuYuAdjYd*Tfd - 3*traceYeAdjYeYeAdjYe*Tfd - tracefuAdjfu*AbsSqr
      (Lambdax)*Tfd - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*Tfd - 2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*Tfd - 3*traceYuAdjYu*AbsSqr(
      Lambdax)*Tfd - AbsSqr(Lambdax)*AbsSqr(Sigmax)*Tfd - 0.4*traceYdAdjYd*Sqr(
      g1)*Tfd + 1.2*traceYeAdjYe*Sqr(g1)*Tfd + tracefdAdjfd*Sqr(g1p)*Tfd - 0.6*
      traceYdAdjYd*Sqr(g1p)*Tfd - 0.2*traceYeAdjYe*Sqr(g1p)*Tfd - 0.25*AbsSqr(
      Lambdax)*Sqr(g1p)*Tfd + 0.27*Sqr(g1)*Sqr(g1p)*Tfd + 1.8*Sqr(g1)*Sqr(g2)*
      Tfd + 1.95*Sqr(g1p)*Sqr(g2)*Tfd + 16*traceYdAdjYd*Sqr(g3)*Tfd + 0.095*
      Power(g1p,4)*Sqr(QS)*Tfd + 0.05*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*Tfd - 3*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)*Tfd - 0.02*fd*(1188*Power(g1,4)*MassB +
      3933*Power(g1p,4)*MassBp + 3300*Power(g2,4)*MassWB + 100*
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
      (g1p) - 20*MassBp*traceYeAdjYe*Sqr(g1p) + 27*MassB*Sqr(g1)*Sqr(g1p) + 27*
      MassBp*Sqr(g1)*Sqr(g1p) + 180*MassB*Sqr(g1)*Sqr(g2) + 180*MassWB*Sqr(g1)*
      Sqr(g2) + 195*MassBp*Sqr(g1p)*Sqr(g2) + 195*MassWB*Sqr(g1p)*Sqr(g2) -
      1600*traceAdjYdTYd*Sqr(g3) + 1600*MassG*traceYdAdjYd*Sqr(g3) + 19*Power(
      g1p,4)*MassBp*Sqr(QS) + 600*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 5*Conj(
      Lambdax)*((20*tracefuAdjfu + 60*traceKappaAdjKappa + 40*
      traceLambda12AdjLambda12 + 60*traceYuAdjYu + 20*AbsSqr(Sigmax) + 5*Sqr(
      g1p) - Sqr(g1p)*Sqr(QS))*TLambdax + Lambdax*(20*traceAdjfuTfu + 60*
      traceAdjKappaTKappa + 40*traceAdjLambda12TLambda12 + 60*traceAdjYuTYu - 5
      *MassBp*Sqr(g1p) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) -
      0.6*(10*traceAdjfdTfd + 30*traceAdjYdTYd + 10*traceAdjYeTYe + 4*MassB*
      Sqr(g1) + MassBp*Sqr(g1p) + 20*MassWB*Sqr(g2) + 10*Conj(Lambdax)*TLambdax
      )*(fd*fd.adjoint()*fd) - 4*tracefdAdjfd*(fd*fd.adjoint()*Tfd) - 12*
      traceYdAdjYd*(fd*fd.adjoint()*Tfd) - 4*traceYeAdjYe*(fd*fd.adjoint()*Tfd)
      - 4*AbsSqr(Lambdax)*(fd*fd.adjoint()*Tfd) + 1.2*Sqr(g1)*(fd*fd.adjoint()
      *Tfd) + 1.8*Sqr(g1p)*(fd*fd.adjoint()*Tfd) + 6*Sqr(g2)*(fd*fd.adjoint()*
      Tfd) - 6*traceAdjKappaTKappa*(fd*Lambda12.conjugate()*(Lambda12)
      .transpose()) - 4*traceAdjLambda12TLambda12*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()) - 0.5*MassBp*Sqr(g1p)*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()) - 0.1*MassBp*Sqr(g1p)*Sqr(QS)*(fd*
      Lambda12.conjugate()*(Lambda12).transpose()) - 4*Conj(Lambdax)*TLambdax*(
      fd*Lambda12.conjugate()*(Lambda12).transpose()) - 2*Conj(Sigmax)*TSigmax*
      (fd*Lambda12.conjugate()*(Lambda12).transpose()) - 6*traceKappaAdjKappa*(
      fd*Lambda12.conjugate()*(TLambda12).transpose()) - 4*
      traceLambda12AdjLambda12*(fd*Lambda12.conjugate()*(TLambda12).transpose()
      ) - 4*AbsSqr(Lambdax)*(fd*Lambda12.conjugate()*(TLambda12).transpose()) -
      2*AbsSqr(Sigmax)*(fd*Lambda12.conjugate()*(TLambda12).transpose()) + 0.5
      *Sqr(g1p)*(fd*Lambda12.conjugate()*(TLambda12).transpose()) + 0.1*Sqr(g1p
      )*Sqr(QS)*(fd*Lambda12.conjugate()*(TLambda12).transpose()) - 4*
      traceAdjfuTfu*(fu*fu.adjoint()*fd) - 12*traceAdjYuTYu*(fu*fu.adjoint()*fd
      ) - 2.4*MassB*Sqr(g1)*(fu*fu.adjoint()*fd) + 2.4*MassBp*Sqr(g1p)*(fu*
      fu.adjoint()*fd) - 12*MassWB*Sqr(g2)*(fu*fu.adjoint()*fd) - 4*Conj(
      Lambdax)*TLambdax*(fu*fu.adjoint()*fd) - 2*tracefuAdjfu*(fu*fu.adjoint()*
      Tfd) - 6*traceYuAdjYu*(fu*fu.adjoint()*Tfd) - 2*AbsSqr(Lambdax)*(fu*
      fu.adjoint()*Tfd) + 1.2*Sqr(g1)*(fu*fu.adjoint()*Tfd) - 1.2*Sqr(g1p)*(fu*
      fu.adjoint()*Tfd) + 6*Sqr(g2)*(fu*fu.adjoint()*Tfd) - 5*tracefdAdjfd*(Tfd
      *fd.adjoint()*fd) - 15*traceYdAdjYd*(Tfd*fd.adjoint()*fd) - 5*
      traceYeAdjYe*(Tfd*fd.adjoint()*fd) - 5*AbsSqr(Lambdax)*(Tfd*fd.adjoint()*
      fd) + 2.4*Sqr(g1)*(Tfd*fd.adjoint()*fd) - 0.9*Sqr(g1p)*(Tfd*fd.adjoint()*
      fd) + 12*Sqr(g2)*(Tfd*fd.adjoint()*fd) - 3*traceKappaAdjKappa*(Tfd*
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*traceLambda12AdjLambda12
      *(Tfd*Lambda12.conjugate()*(Lambda12).transpose()) - 2*AbsSqr(Lambdax)*(
      Tfd*Lambda12.conjugate()*(Lambda12).transpose()) - AbsSqr(Sigmax)*(Tfd*
      Lambda12.conjugate()*(Lambda12).transpose()) + 0.25*Sqr(g1p)*(Tfd*
      Lambda12.conjugate()*(Lambda12).transpose()) + 0.05*Sqr(g1p)*Sqr(QS)*(Tfd
      *Lambda12.conjugate()*(Lambda12).transpose()) - 4*tracefuAdjfu*(Tfu*
      fu.adjoint()*fd) - 12*traceYuAdjYu*(Tfu*fu.adjoint()*fd) - 4*AbsSqr(
      Lambdax)*(Tfu*fu.adjoint()*fd) + 2.4*Sqr(g1)*(Tfu*fu.adjoint()*fd) - 2.4*
      Sqr(g1p)*(Tfu*fu.adjoint()*fd) + 12*Sqr(g2)*(Tfu*fu.adjoint()*fd) - 6*(fd
      *fd.adjoint()*fd*fd.adjoint()*Tfd) - 4*(fd*fd.adjoint()*fu*fu.adjoint()*
      Tfd) - 8*(fd*fd.adjoint()*Tfd*fd.adjoint()*fd) - 4*(fd*fd.adjoint()*Tfu*
      fu.adjoint()*fd) - 2*(fd*Lambda12.conjugate()*fu.transpose()*fu.conjugate
      ()*(TLambda12).transpose()) - 2*(fd*Lambda12.conjugate()*hE.transpose()*
      hE.conjugate()*(TLambda12).transpose()) - 2*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()*fd.adjoint()*Tfd) - 2*(fd*Lambda12.conjugate()*(
      Lambda12).transpose()*Lambda12.conjugate()*(TLambda12).transpose()) - 2*(
      fd*Lambda12.conjugate()*(Tfu).transpose()*fu.conjugate()*(Lambda12)
      .transpose()) - 2*(fd*Lambda12.conjugate()*(ThE).transpose()*hE.conjugate
      ()*(Lambda12).transpose()) - 4*(fd*Lambda12.conjugate()*(TLambda12)
      .transpose()*fd.adjoint()*fd) - 2*(fd*Lambda12.conjugate()*(TLambda12)
      .transpose()*Lambda12.conjugate()*(Lambda12).transpose()) - 2*(fu*
      fu.adjoint()*fu*fu.adjoint()*Tfd) - 4*(fu*fu.adjoint()*Tfu*fu.adjoint()*
      fd) - 2*(fu*hE.adjoint()*hE*fu.adjoint()*Tfd) - 4*(fu*hE.adjoint()*ThE*
      fu.adjoint()*fd) - 2*(fu*(Lambda12).adjoint()*Lambda12*fu.adjoint()*Tfd)
      - 4*(fu*(Lambda12).adjoint()*TLambda12*fu.adjoint()*fd) - 6*(Tfd*
      fd.adjoint()*fd*fd.adjoint()*fd) - 2*(Tfd*fd.adjoint()*fu*fu.adjoint()*fd
      ) - Tfd*Lambda12.conjugate()*fu.transpose()*fu.conjugate()*(Lambda12)
      .transpose() - Tfd*Lambda12.conjugate()*hE.transpose()*hE.conjugate()*(
      Lambda12).transpose() - 4*(Tfd*Lambda12.conjugate()*(Lambda12).transpose(
      )*fd.adjoint()*fd) - Tfd*Lambda12.conjugate()*(Lambda12).transpose()*
      Lambda12.conjugate()*(Lambda12).transpose() - 4*(Tfu*fu.adjoint()*fu*
      fu.adjoint()*fd) - 4*(Tfu*hE.adjoint()*hE*fu.adjoint()*fd) - 4*(Tfu*(
      Lambda12).adjoint()*Lambda12*fu.adjoint()*fd))).real();


   return beta_Tfd;
}

} // namespace flexiblesusy
