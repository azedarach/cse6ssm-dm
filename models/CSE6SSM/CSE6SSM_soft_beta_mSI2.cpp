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

// File generated at Wed 3 Jun 2015 23:43:36

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mSI2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mSI2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mSI2;

   beta_mSI2 = (oneOver16PiSqr*(2*(2*mHd2*(fd.conjugate()*fd.transpose())
      + 2*mHu2*(fu.conjugate()*fu.transpose()) + 2*(Tfd.conjugate()*(Tfd)
      .transpose()) + 2*(Tfu.conjugate()*(Tfu).transpose()) + mSI2*fd.conjugate
      ()*fd.transpose() + mSI2*fu.conjugate()*fu.transpose() + 2*(fd.conjugate(
      )*mH2I2.conjugate()*fd.transpose()) + fd.conjugate()*fd.transpose()*mSI2
      + 2*(fu.conjugate()*mH1I2.conjugate()*fu.transpose()) + fu.conjugate()*
      fu.transpose()*mSI2) + 1.5811388300841898*g1p*Tr14*UNITMATRIX(3) - 5*
      AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(3))).real();


   return beta_mSI2;
}

/**
 * Calculates the two-loop beta function of mSI2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mSI2_two_loop(const Soft_traces& soft_traces) const
{
   const double tracefdAdjfd = TRACE_STRUCT.tracefdAdjfd;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTfdTpTfd = TRACE_STRUCT.traceconjTfdTpTfd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracefdmH2I2Adjfd = TRACE_STRUCT.tracefdmH2I2Adjfd;
   const double tracefdAdjfdconjmSI2 = TRACE_STRUCT.tracefdAdjfdconjmSI2;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTfdTpfd = TRACE_STRUCT.traceconjTfdTpfd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double tracefuAdjfu = TRACE_STRUCT.tracefuAdjfu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTfuTpTfu = TRACE_STRUCT.traceconjTfuTpTfu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracefumH1I2Adjfu = TRACE_STRUCT.tracefumH1I2Adjfu;
   const double tracefuAdjfuconjmSI2 = TRACE_STRUCT.tracefuAdjfuconjmSI2;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceconjTfuTpfu = TRACE_STRUCT.traceconjTfuTpfu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mSI2;

   beta_mSI2 = (twoLoop*(-4*traceconjTfdTpTfd*(fd.conjugate()*
      fd.transpose()) - 12*traceconjTYdTpTYd*(fd.conjugate()*fd.transpose()) -
      4*traceconjTYeTpTYe*(fd.conjugate()*fd.transpose()) - 8*mHd2*tracefdAdjfd
      *(fd.conjugate()*fd.transpose()) - 4*tracefdAdjfdconjmSI2*(fd.conjugate()
      *fd.transpose()) - 4*tracefdmH2I2Adjfd*(fd.conjugate()*fd.transpose()) -
      12*tracemd2YdAdjYd*(fd.conjugate()*fd.transpose()) - 4*traceme2YeAdjYe*(
      fd.conjugate()*fd.transpose()) - 4*traceml2AdjYeYe*(fd.conjugate()*
      fd.transpose()) - 12*tracemq2AdjYdYd*(fd.conjugate()*fd.transpose()) - 24
      *mHd2*traceYdAdjYd*(fd.conjugate()*fd.transpose()) - 8*mHd2*traceYeAdjYe*
      (fd.conjugate()*fd.transpose()) - 8*mHd2*AbsSqr(Lambdax)*(fd.conjugate()*
      fd.transpose()) - 4*mHu2*AbsSqr(Lambdax)*(fd.conjugate()*fd.transpose())
      - 4*ms2*AbsSqr(Lambdax)*(fd.conjugate()*fd.transpose()) - 4*AbsSqr(
      TLambdax)*(fd.conjugate()*fd.transpose()) + 2.4*mHd2*Sqr(g1)*(
      fd.conjugate()*fd.transpose()) + 4.8*AbsSqr(MassB)*Sqr(g1)*(fd.conjugate(
      )*fd.transpose()) - 2.4*mHd2*Sqr(g1p)*(fd.conjugate()*fd.transpose()) +
      12*mHd2*Sqr(g2)*(fd.conjugate()*fd.transpose()) + 24*AbsSqr(MassWB)*Sqr(
      g2)*(fd.conjugate()*fd.transpose()) - 4*traceconjTfdTpfd*(fd.conjugate()*
      (Tfd).transpose()) - 12*traceconjTYdTpYd*(fd.conjugate()*(Tfd).transpose(
      )) - 4*traceconjTYeTpYe*(fd.conjugate()*(Tfd).transpose()) - 4*Conj(
      TLambdax)*Lambdax*(fd.conjugate()*(Tfd).transpose()) - 2.4*Conj(MassB)*
      Sqr(g1)*(fd.conjugate()*(Tfd).transpose()) - 12*Conj(MassWB)*Sqr(g2)*(
      fd.conjugate()*(Tfd).transpose()) - 4*traceconjTfuTpTfu*(fu.conjugate()*
      fu.transpose()) - 12*traceconjTYuTpTYu*(fu.conjugate()*fu.transpose()) -
      8*mHu2*tracefuAdjfu*(fu.conjugate()*fu.transpose()) - 4*
      tracefuAdjfuconjmSI2*(fu.conjugate()*fu.transpose()) - 4*
      tracefumH1I2Adjfu*(fu.conjugate()*fu.transpose()) - 12*tracemq2AdjYuYu*(
      fu.conjugate()*fu.transpose()) - 12*tracemu2YuAdjYu*(fu.conjugate()*
      fu.transpose()) - 24*mHu2*traceYuAdjYu*(fu.conjugate()*fu.transpose()) -
      4*mHd2*AbsSqr(Lambdax)*(fu.conjugate()*fu.transpose()) - 8*mHu2*AbsSqr(
      Lambdax)*(fu.conjugate()*fu.transpose()) - 4*ms2*AbsSqr(Lambdax)*(
      fu.conjugate()*fu.transpose()) - 4*AbsSqr(TLambdax)*(fu.conjugate()*
      fu.transpose()) + 2.4*mHu2*Sqr(g1)*(fu.conjugate()*fu.transpose()) + 4.8*
      AbsSqr(MassB)*Sqr(g1)*(fu.conjugate()*fu.transpose()) - 2.4*mHu2*Sqr(g1p)
      *(fu.conjugate()*fu.transpose()) + 12*mHu2*Sqr(g2)*(fu.conjugate()*
      fu.transpose()) + 24*AbsSqr(MassWB)*Sqr(g2)*(fu.conjugate()*fu.transpose(
      )) - 4*traceconjTfuTpfu*(fu.conjugate()*(Tfu).transpose()) - 12*
      traceconjTYuTpYu*(fu.conjugate()*(Tfu).transpose()) - 4*Conj(TLambdax)*
      Lambdax*(fu.conjugate()*(Tfu).transpose()) - 2.4*Conj(MassB)*Sqr(g1)*(
      fu.conjugate()*(Tfu).transpose()) - 12*Conj(MassWB)*Sqr(g2)*(fu.conjugate
      ()*(Tfu).transpose()) - 4*traceAdjfdTfd*(Tfd.conjugate()*fd.transpose())
      - 12*traceAdjYdTYd*(Tfd.conjugate()*fd.transpose()) - 4*traceAdjYeTYe*(
      Tfd.conjugate()*fd.transpose()) - 2.4*MassB*Sqr(g1)*(Tfd.conjugate()*
      fd.transpose()) + 2.4*MassBp*Sqr(g1p)*(Tfd.conjugate()*fd.transpose()) -
      12*MassWB*Sqr(g2)*(Tfd.conjugate()*fd.transpose()) - 4*Conj(Lambdax)*
      TLambdax*(Tfd.conjugate()*fd.transpose()) - 4*tracefdAdjfd*(Tfd.conjugate
      ()*(Tfd).transpose()) - 12*traceYdAdjYd*(Tfd.conjugate()*(Tfd).transpose(
      )) - 4*traceYeAdjYe*(Tfd.conjugate()*(Tfd).transpose()) - 4*AbsSqr(
      Lambdax)*(Tfd.conjugate()*(Tfd).transpose()) + 2.4*Sqr(g1)*(Tfd.conjugate
      ()*(Tfd).transpose()) - 2.4*Sqr(g1p)*(Tfd.conjugate()*(Tfd).transpose())
      + 12*Sqr(g2)*(Tfd.conjugate()*(Tfd).transpose()) - 4*traceAdjfuTfu*(
      Tfu.conjugate()*fu.transpose()) - 12*traceAdjYuTYu*(Tfu.conjugate()*
      fu.transpose()) - 2.4*MassB*Sqr(g1)*(Tfu.conjugate()*fu.transpose()) +
      2.4*MassBp*Sqr(g1p)*(Tfu.conjugate()*fu.transpose()) - 12*MassWB*Sqr(g2)*
      (Tfu.conjugate()*fu.transpose()) - 4*Conj(Lambdax)*TLambdax*(
      Tfu.conjugate()*fu.transpose()) - 4*tracefuAdjfu*(Tfu.conjugate()*(Tfu)
      .transpose()) - 12*traceYuAdjYu*(Tfu.conjugate()*(Tfu).transpose()) - 4*
      AbsSqr(Lambdax)*(Tfu.conjugate()*(Tfu).transpose()) + 2.4*Sqr(g1)*(
      Tfu.conjugate()*(Tfu).transpose()) - 2.4*Sqr(g1p)*(Tfu.conjugate()*(Tfu)
      .transpose()) + 12*Sqr(g2)*(Tfu.conjugate()*(Tfu).transpose()) - 2*
      tracefdAdjfd*(mSI2*fd.conjugate()*fd.transpose()) - 6*traceYdAdjYd*(mSI2*
      fd.conjugate()*fd.transpose()) - 2*traceYeAdjYe*(mSI2*fd.conjugate()*
      fd.transpose()) - 2*AbsSqr(Lambdax)*(mSI2*fd.conjugate()*fd.transpose())
      + 1.2*Sqr(g1)*(mSI2*fd.conjugate()*fd.transpose()) - 1.2*Sqr(g1p)*(mSI2*
      fd.conjugate()*fd.transpose()) + 6*Sqr(g2)*(mSI2*fd.conjugate()*
      fd.transpose()) - 2*tracefuAdjfu*(mSI2*fu.conjugate()*fu.transpose()) - 6
      *traceYuAdjYu*(mSI2*fu.conjugate()*fu.transpose()) - 2*AbsSqr(Lambdax)*(
      mSI2*fu.conjugate()*fu.transpose()) + 1.2*Sqr(g1)*(mSI2*fu.conjugate()*
      fu.transpose()) - 1.2*Sqr(g1p)*(mSI2*fu.conjugate()*fu.transpose()) + 6*
      Sqr(g2)*(mSI2*fu.conjugate()*fu.transpose()) - 4*tracefdAdjfd*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 12*traceYdAdjYd*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 4*traceYeAdjYe*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 4*AbsSqr(Lambdax)*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) + 2.4*Sqr(g1)*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 2.4*Sqr(g1p)*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) + 12*Sqr(g2)*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 2*tracefdAdjfd*(
      fd.conjugate()*fd.transpose()*mSI2) - 6*traceYdAdjYd*(fd.conjugate()*
      fd.transpose()*mSI2) - 2*traceYeAdjYe*(fd.conjugate()*fd.transpose()*mSI2
      ) - 2*AbsSqr(Lambdax)*(fd.conjugate()*fd.transpose()*mSI2) + 1.2*Sqr(g1)*
      (fd.conjugate()*fd.transpose()*mSI2) - 1.2*Sqr(g1p)*(fd.conjugate()*
      fd.transpose()*mSI2) + 6*Sqr(g2)*(fd.conjugate()*fd.transpose()*mSI2) - 4
      *tracefuAdjfu*(fu.conjugate()*mH1I2.conjugate()*fu.transpose()) - 12*
      traceYuAdjYu*(fu.conjugate()*mH1I2.conjugate()*fu.transpose()) - 4*AbsSqr
      (Lambdax)*(fu.conjugate()*mH1I2.conjugate()*fu.transpose()) + 2.4*Sqr(g1)
      *(fu.conjugate()*mH1I2.conjugate()*fu.transpose()) - 2.4*Sqr(g1p)*(
      fu.conjugate()*mH1I2.conjugate()*fu.transpose()) + 12*Sqr(g2)*(
      fu.conjugate()*mH1I2.conjugate()*fu.transpose()) - 2*tracefuAdjfu*(
      fu.conjugate()*fu.transpose()*mSI2) - 6*traceYuAdjYu*(fu.conjugate()*
      fu.transpose()*mSI2) - 2*AbsSqr(Lambdax)*(fu.conjugate()*fu.transpose()*
      mSI2) + 1.2*Sqr(g1)*(fu.conjugate()*fu.transpose()*mSI2) - 1.2*Sqr(g1p)*(
      fu.conjugate()*fu.transpose()*mSI2) + 6*Sqr(g2)*(fu.conjugate()*
      fu.transpose()*mSI2) - 4*mHd2*(fd.conjugate()*Lambda12*(Lambda12).adjoint
      ()*fd.transpose()) - 4*ms2*(fd.conjugate()*Lambda12*(Lambda12).adjoint()*
      fd.transpose()) - 4*(fd.conjugate()*Lambda12*(TLambda12).adjoint()*(Tfd)
      .transpose()) - 4*(fd.conjugate()*TLambda12*(TLambda12).adjoint()*
      fd.transpose()) - 8*mHd2*(fd.conjugate()*fd.transpose()*fd.conjugate()*
      fd.transpose()) - 4*(fd.conjugate()*fd.transpose()*Tfd.conjugate()*(Tfd)
      .transpose()) - 4*(fd.conjugate()*(Tfd).transpose()*Tfd.conjugate()*
      fd.transpose()) - 8*mHu2*(fu.conjugate()*fu.transpose()*fu.conjugate()*
      fu.transpose()) - 4*(fu.conjugate()*fu.transpose()*Tfu.conjugate()*(Tfu)
      .transpose()) - 4*mHp2*(fu.conjugate()*hE.transpose()*hE.conjugate()*
      fu.transpose()) - 4*mHu2*(fu.conjugate()*hE.transpose()*hE.conjugate()*
      fu.transpose()) - 4*(fu.conjugate()*hE.transpose()*ThE.conjugate()*(Tfu)
      .transpose()) - 4*mHu2*(fu.conjugate()*(Lambda12).transpose()*
      Lambda12.conjugate()*fu.transpose()) - 4*ms2*(fu.conjugate()*(Lambda12)
      .transpose()*Lambda12.conjugate()*fu.transpose()) - 4*(fu.conjugate()*(
      Lambda12).transpose()*TLambda12.conjugate()*(Tfu).transpose()) - 4*(
      fu.conjugate()*(Tfu).transpose()*Tfu.conjugate()*fu.transpose()) - 4*(
      fu.conjugate()*(ThE).transpose()*ThE.conjugate()*fu.transpose()) - 4*(
      fu.conjugate()*(TLambda12).transpose()*TLambda12.conjugate()*fu.transpose
      ()) - 4*(Tfd.conjugate()*Lambda12*(Lambda12).adjoint()*(Tfd).transpose())
      - 4*(Tfd.conjugate()*TLambda12*(Lambda12).adjoint()*fd.transpose()) - 4*
      (Tfd.conjugate()*fd.transpose()*fd.conjugate()*(Tfd).transpose()) - 4*(
      Tfd.conjugate()*(Tfd).transpose()*fd.conjugate()*fd.transpose()) - 4*(
      Tfu.conjugate()*fu.transpose()*fu.conjugate()*(Tfu).transpose()) - 4*(
      Tfu.conjugate()*hE.transpose()*hE.conjugate()*(Tfu).transpose()) - 4*(
      Tfu.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(Tfu)
      .transpose()) - 4*(Tfu.conjugate()*(Tfu).transpose()*fu.conjugate()*
      fu.transpose()) - 4*(Tfu.conjugate()*(ThE).transpose()*hE.conjugate()*
      fu.transpose()) - 4*(Tfu.conjugate()*(TLambda12).transpose()*
      Lambda12.conjugate()*fu.transpose()) - 2*(mSI2*fd.conjugate()*Lambda12*(
      Lambda12).adjoint()*fd.transpose()) - 2*(mSI2*fd.conjugate()*fd.transpose
      ()*fd.conjugate()*fd.transpose()) - 2*(mSI2*fu.conjugate()*fu.transpose()
      *fu.conjugate()*fu.transpose()) - 2*(mSI2*fu.conjugate()*hE.transpose()*
      hE.conjugate()*fu.transpose()) - 2*(mSI2*fu.conjugate()*(Lambda12)
      .transpose()*Lambda12.conjugate()*fu.transpose()) - 4*(fd.conjugate()*
      mH2I2.conjugate()*Lambda12*(Lambda12).adjoint()*fd.transpose()) - 4*(
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()*fd.conjugate()*
      fd.transpose()) - 4*(fd.conjugate()*Lambda12*mH1I2*(Lambda12).adjoint()*
      fd.transpose()) - 4*(fd.conjugate()*Lambda12*(Lambda12).adjoint()*
      mH2I2.conjugate()*fd.transpose()) - 2*(fd.conjugate()*Lambda12*(Lambda12)
      .adjoint()*fd.transpose()*mSI2) - 4*(fd.conjugate()*fd.transpose()*mSI2*
      fd.conjugate()*fd.transpose()) - 4*(fd.conjugate()*fd.transpose()*
      fd.conjugate()*mH2I2.conjugate()*fd.transpose()) - 2*(fd.conjugate()*
      fd.transpose()*fd.conjugate()*fd.transpose()*mSI2) - 4*(fu.conjugate()*
      mH1I2.conjugate()*fu.transpose()*fu.conjugate()*fu.transpose()) - 4*(
      fu.conjugate()*mH1I2.conjugate()*hE.transpose()*hE.conjugate()*
      fu.transpose()) - 4*(fu.conjugate()*mH1I2.conjugate()*(Lambda12)
      .transpose()*Lambda12.conjugate()*fu.transpose()) - 4*(fu.conjugate()*
      fu.transpose()*mSI2*fu.conjugate()*fu.transpose()) - 4*(fu.conjugate()*
      fu.transpose()*fu.conjugate()*mH1I2.conjugate()*fu.transpose()) - 2*(
      fu.conjugate()*fu.transpose()*fu.conjugate()*fu.transpose()*mSI2) - 4*(
      fu.conjugate()*hE.transpose()*hE.conjugate()*mH1I2.conjugate()*
      fu.transpose()) - 2*(fu.conjugate()*hE.transpose()*hE.conjugate()*
      fu.transpose()*mSI2) - 4*(fu.conjugate()*hE.transpose()*me2.conjugate()*
      hE.conjugate()*fu.transpose()) - 4*(fu.conjugate()*(Lambda12).transpose()
      *mH2I2*Lambda12.conjugate()*fu.transpose()) - 4*(fu.conjugate()*(Lambda12
      ).transpose()*Lambda12.conjugate()*mH1I2.conjugate()*fu.transpose()) - 2*
      (fu.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*fu.transpose(
      )*mSI2) + g1p*(5*g1p*Tr2U144 + 6.324555320336759*Tr34)*UNITMATRIX(3) +
      0.15*Conj(MassBp)*Sqr(g1p)*(16*(-2*MassBp*(fd.conjugate()*fd.transpose())
      + fd.conjugate()*(Tfd).transpose() - 2*MassBp*(fu.conjugate()*
      fu.transpose()) + fu.conjugate()*(Tfu).transpose()) + 5*MassBp*Sqr(g1p)*(
      213 + Sqr(QS))*UNITMATRIX(3)))).real();


   return beta_mSI2;
}

} // namespace flexiblesusy
