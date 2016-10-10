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

// File generated at Wed 3 Jun 2015 23:43:27

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu)
      .adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*
      Yu.adjoint()*mu2) - 1.0327955589886444*g1*Tr11*UNITMATRIX(3) +
      0.31622776601683794*g1p*Tr14*UNITMATRIX(3) - 2.1333333333333333*AbsSqr(
      MassB)*Sqr(g1)*UNITMATRIX(3) - 0.2*AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(3)
      - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTfuTpTfu = TRACE_STRUCT.traceconjTfuTpTfu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracefumH1I2Adjfu = TRACE_STRUCT.tracefumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2 = TRACE_STRUCT.tracefuAdjfuconjmSI2;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTfuTpfu = TRACE_STRUCT.traceconjTfuTpfu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (twoLoop*(-4*traceconjTfuTpTfu*(Yu*Yu.adjoint()) - 12*
      traceconjTYuTpTYu*(Yu*Yu.adjoint()) - 8*mHu2*tracefuAdjfu*(Yu*Yu.adjoint(
      )) - 4*tracefuAdjfuconjmSI2*(Yu*Yu.adjoint()) - 4*tracefumH1I2Adjfu*(Yu*
      Yu.adjoint()) - 12*tracemq2AdjYuYu*(Yu*Yu.adjoint()) - 12*tracemu2YuAdjYu
      *(Yu*Yu.adjoint()) - 24*mHu2*traceYuAdjYu*(Yu*Yu.adjoint()) - 4*mHd2*
      AbsSqr(Lambdax)*(Yu*Yu.adjoint()) - 8*mHu2*AbsSqr(Lambdax)*(Yu*Yu.adjoint
      ()) - 4*ms2*AbsSqr(Lambdax)*(Yu*Yu.adjoint()) - 4*AbsSqr(TLambdax)*(Yu*
      Yu.adjoint()) - 0.8*mHu2*Sqr(g1)*(Yu*Yu.adjoint()) + 0.8*mHu2*Sqr(g1p)*(
      Yu*Yu.adjoint()) + 12*mHu2*Sqr(g2)*(Yu*Yu.adjoint()) + 24*AbsSqr(MassWB)*
      Sqr(g2)*(Yu*Yu.adjoint()) - 4*traceAdjfuTfu*(Yu*(TYu).adjoint()) - 12*
      traceAdjYuTYu*(Yu*(TYu).adjoint()) + 0.8*MassB*Sqr(g1)*(Yu*(TYu).adjoint(
      )) - 0.8*MassBp*Sqr(g1p)*(Yu*(TYu).adjoint()) - 12*MassWB*Sqr(g2)*(Yu*(
      TYu).adjoint()) - 4*Conj(Lambdax)*TLambdax*(Yu*(TYu).adjoint()) - 4*
      traceconjTfuTpfu*(TYu*Yu.adjoint()) - 12*traceconjTYuTpYu*(TYu*Yu.adjoint
      ()) - 4*Conj(TLambdax)*Lambdax*(TYu*Yu.adjoint()) - 12*Conj(MassWB)*Sqr(
      g2)*(TYu*Yu.adjoint()) - 4*tracefuAdjfu*(TYu*(TYu).adjoint()) - 12*
      traceYuAdjYu*(TYu*(TYu).adjoint()) - 4*AbsSqr(Lambdax)*(TYu*(TYu).adjoint
      ()) - 0.8*Sqr(g1)*(TYu*(TYu).adjoint()) + 0.8*Sqr(g1p)*(TYu*(TYu).adjoint
      ()) + 12*Sqr(g2)*(TYu*(TYu).adjoint()) - 2*tracefuAdjfu*(mu2*Yu*
      Yu.adjoint()) - 6*traceYuAdjYu*(mu2*Yu*Yu.adjoint()) - 2*AbsSqr(Lambdax)*
      (mu2*Yu*Yu.adjoint()) - 0.4*Sqr(g1)*(mu2*Yu*Yu.adjoint()) + 0.4*Sqr(g1p)*
      (mu2*Yu*Yu.adjoint()) + 6*Sqr(g2)*(mu2*Yu*Yu.adjoint()) - 4*tracefuAdjfu*
      (Yu*mq2*Yu.adjoint()) - 12*traceYuAdjYu*(Yu*mq2*Yu.adjoint()) - 4*AbsSqr(
      Lambdax)*(Yu*mq2*Yu.adjoint()) - 0.8*Sqr(g1)*(Yu*mq2*Yu.adjoint()) + 0.8*
      Sqr(g1p)*(Yu*mq2*Yu.adjoint()) + 12*Sqr(g2)*(Yu*mq2*Yu.adjoint()) - 2*
      tracefuAdjfu*(Yu*Yu.adjoint()*mu2) - 6*traceYuAdjYu*(Yu*Yu.adjoint()*mu2)
      - 2*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*mu2) - 0.4*Sqr(g1)*(Yu*Yu.adjoint()
      *mu2) + 0.4*Sqr(g1p)*(Yu*Yu.adjoint()*mu2) + 6*Sqr(g2)*(Yu*Yu.adjoint()*
      mu2) - 4*mHd2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 4*mHu2*(Yu*Yd.adjoint()
      *Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*TYd*(TYu).adjoint()) - 8*mHu2*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yu.adjoint()*TYu*(TYu).adjoint()) -
      4*(Yu*(TYd).adjoint()*TYd*Yu.adjoint()) - 4*(Yu*(TYu).adjoint()*TYu*
      Yu.adjoint()) - 4*mHp2*(Yu*gD.conjugate()*gD.transpose()*Yu.adjoint()) -
      4*mHu2*(Yu*gD.conjugate()*gD.transpose()*Yu.adjoint()) - 4*(Yu*
      gD.conjugate()*(TgD).transpose()*(TYu).adjoint()) - 4*(Yu*TgD.conjugate()
      *(TgD).transpose()*Yu.adjoint()) - 4*(TYu*Yd.adjoint()*Yd*(TYu).adjoint()
      ) - 4*(TYu*Yu.adjoint()*Yu*(TYu).adjoint()) - 4*(TYu*(TYd).adjoint()*Yd*
      Yu.adjoint()) - 4*(TYu*(TYu).adjoint()*Yu*Yu.adjoint()) - 4*(TYu*
      gD.conjugate()*gD.transpose()*(TYu).adjoint()) - 4*(TYu*TgD.conjugate()*
      gD.transpose()*Yu.adjoint()) - 2*(mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) -
      2*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 2*(mu2*Yu*gD.conjugate()*
      gD.transpose()*Yu.adjoint()) - 4*(Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) -
      4*(Yu*mq2*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*mq2*gD.conjugate()*
      gD.transpose()*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) -
      4*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*(Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*mu2) - 4*(Yu*Yu.adjoint()*mu2*Yu*Yu.adjoint()) - 4*(Yu*
      Yu.adjoint()*Yu*mq2*Yu.adjoint()) - 2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*
      mu2) - 4*(Yu*gD.conjugate()*mDxbar2*gD.transpose()*Yu.adjoint()) - 4*(Yu*
      gD.conjugate()*gD.transpose()*mq2*Yu.adjoint()) - 2*(Yu*gD.conjugate()*
      gD.transpose()*Yu.adjoint()*mu2) + 10.666666666666666*Power(g3,4)*Tr23*
      UNITMATRIX(3) - 0.6531972647421809*g1*g1p*Tr2U114*UNITMATRIX(3) -
      0.6531972647421809*g1*g1p*Tr2U141*UNITMATRIX(3) - 4.131182235954578*g1*
      Tr31*UNITMATRIX(3) + 1.2649110640673518*g1p*Tr34*UNITMATRIX(3) +
      53.333333333333336*Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) +
      2.1333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 0.2*Tr2U144*Sqr(g1p)*
      UNITMATRIX(3) + 11.377777777777778*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 5.688888888888889*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 1.0666666666666667*AbsSqr(MassG)*Sqr(g1p)*Sqr(g3)*
      UNITMATRIX(3) + 0.5333333333333333*MassBp*Conj(MassG)*Sqr(g1p)*Sqr(g3)*
      UNITMATRIX(3) + 0.017777777777777778*Conj(MassB)*Sqr(g1)*(45*(-2*MassB*(
      Yu*Yu.adjoint()) + TYu*Yu.adjoint()) + 8*(456*MassB*Sqr(g1) + 3*(2*MassB
      + MassBp)*Sqr(g1p) + 40*(2*MassB + MassG)*Sqr(g3))*UNITMATRIX(3)) +
      0.0033333333333333335*Conj(MassBp)*Sqr(g1p)*(240*(2*MassBp*(Yu*Yu.adjoint
      ()) - TYu*Yu.adjoint()) + (128*(MassB + 2*MassBp)*Sqr(g1) + 160*(2*MassBp
      + MassG)*Sqr(g3) + 9*MassBp*Sqr(g1p)*(189 + Sqr(QS)))*UNITMATRIX(3))))
      .real();


   return beta_mu2;
}

} // namespace flexiblesusy
