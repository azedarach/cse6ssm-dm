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

// File generated at Wed 3 Jun 2015 23:43:35

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mH2I2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_mH2I2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_mH2I2;

   beta_mH2I2 = (oneOver16PiSqr*(2*mHd2*(fd.adjoint()*fd) + 2*((Tfd)
      .adjoint()*Tfd) + 2*ms2*(Lambda12.conjugate()*(Lambda12).transpose()) + 2
      *(TLambda12.conjugate()*(TLambda12).transpose()) + mH2I2*fd.adjoint()*fd
      + mH2I2*Lambda12.conjugate()*(Lambda12).transpose() + fd.adjoint()*fd*
      mH2I2 + 2*(fd.adjoint()*mSI2.conjugate()*fd) + 2*(Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) + Lambda12.conjugate()*(
      Lambda12).transpose()*mH2I2 + 0.7745966692414834*g1*Tr11*UNITMATRIX(2) -
      0.6324555320336759*g1p*Tr14*UNITMATRIX(2) - 1.2*AbsSqr(MassB)*Sqr(g1)*
      UNITMATRIX(2) - 0.8*AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(2) - 6*AbsSqr(
      MassWB)*Sqr(g2)*UNITMATRIX(2))).real();


   return beta_mH2I2;
}

/**
 * Calculates the two-loop beta function of mH2I2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_mH2I2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,2,2> beta_mH2I2;

   beta_mH2I2 = (twoLoop*(-2*traceconjTfdTpTfd*(fd.adjoint()*fd) - 6*
      traceconjTYdTpTYd*(fd.adjoint()*fd) - 2*traceconjTYeTpTYe*(fd.adjoint()*
      fd) - 4*mHd2*tracefdAdjfd*(fd.adjoint()*fd) - 2*tracefdAdjfdconjmSI2*(
      fd.adjoint()*fd) - 2*tracefdmH2I2Adjfd*(fd.adjoint()*fd) - 6*
      tracemd2YdAdjYd*(fd.adjoint()*fd) - 2*traceme2YeAdjYe*(fd.adjoint()*fd) -
      2*traceml2AdjYeYe*(fd.adjoint()*fd) - 6*tracemq2AdjYdYd*(fd.adjoint()*fd
      ) - 12*mHd2*traceYdAdjYd*(fd.adjoint()*fd) - 4*mHd2*traceYeAdjYe*(
      fd.adjoint()*fd) - 4*mHd2*AbsSqr(Lambdax)*(fd.adjoint()*fd) - 2*mHu2*
      AbsSqr(Lambdax)*(fd.adjoint()*fd) - 2*ms2*AbsSqr(Lambdax)*(fd.adjoint()*
      fd) - 2*AbsSqr(TLambdax)*(fd.adjoint()*fd) + 3*mHd2*Sqr(g1p)*(fd.adjoint(
      )*fd) - 2*traceconjTfdTpfd*(fd.adjoint()*Tfd) - 6*traceconjTYdTpYd*(
      fd.adjoint()*Tfd) - 2*traceconjTYeTpYe*(fd.adjoint()*Tfd) - 2*Conj(
      TLambdax)*Lambdax*(fd.adjoint()*Tfd) - 2*traceAdjfdTfd*((Tfd).adjoint()*
      fd) - 6*traceAdjYdTYd*((Tfd).adjoint()*fd) - 2*traceAdjYeTYe*((Tfd)
      .adjoint()*fd) - 3*MassBp*Sqr(g1p)*((Tfd).adjoint()*fd) - 2*Conj(Lambdax)
      *TLambdax*((Tfd).adjoint()*fd) - 2*tracefdAdjfd*((Tfd).adjoint()*Tfd) - 6
      *traceYdAdjYd*((Tfd).adjoint()*Tfd) - 2*traceYeAdjYe*((Tfd).adjoint()*Tfd
      ) - 2*AbsSqr(Lambdax)*((Tfd).adjoint()*Tfd) + 3*Sqr(g1p)*((Tfd).adjoint()
      *Tfd) - 6*traceconjTKappaTpTKappa*(Lambda12.conjugate()*(Lambda12)
      .transpose()) - 4*traceconjTLambda12TpTLambda12*(Lambda12.conjugate()*(
      Lambda12).transpose()) - 12*ms2*traceKappaAdjKappa*(Lambda12.conjugate()*
      (Lambda12).transpose()) - 6*traceKappaAdjKappaconjmDx2*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 6*
      traceKappaconjmDxbar2AdjKappa*(Lambda12.conjugate()*(Lambda12).transpose(
      )) - 8*ms2*traceLambda12AdjLambda12*(Lambda12.conjugate()*(Lambda12)
      .transpose()) - 4*traceLambda12AdjLambda12conjmH2I2*(Lambda12.conjugate()
      *(Lambda12).transpose()) - 4*tracemH1I2AdjLambda12Lambda12*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 4*mHd2*AbsSqr(Lambdax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 4*mHu2*AbsSqr(Lambdax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 8*ms2*AbsSqr(Lambdax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*mphi2*AbsSqr(Sigmax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 4*ms2*AbsSqr(Sigmax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*msbar2*AbsSqr(Sigmax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 4*AbsSqr(TLambdax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*AbsSqr(TSigmax)*(
      Lambda12.conjugate()*(Lambda12).transpose()) + 0.5*ms2*Sqr(g1p)*(
      Lambda12.conjugate()*(Lambda12).transpose()) + 0.1*ms2*Sqr(g1p)*Sqr(QS)*(
      Lambda12.conjugate()*(Lambda12).transpose()) - 6*traceconjTKappaTpKappa*(
      Lambda12.conjugate()*(TLambda12).transpose()) - 4*
      traceconjTLambda12TpLambda12*(Lambda12.conjugate()*(TLambda12).transpose(
      )) - 4*Conj(TLambdax)*Lambdax*(Lambda12.conjugate()*(TLambda12).transpose
      ()) - 2*Conj(TSigmax)*Sigmax*(Lambda12.conjugate()*(TLambda12).transpose(
      )) - 6*traceAdjKappaTKappa*(TLambda12.conjugate()*(Lambda12).transpose())
      - 4*traceAdjLambda12TLambda12*(TLambda12.conjugate()*(Lambda12)
      .transpose()) - 0.5*MassBp*Sqr(g1p)*(TLambda12.conjugate()*(Lambda12)
      .transpose()) - 0.1*MassBp*Sqr(g1p)*Sqr(QS)*(TLambda12.conjugate()*(
      Lambda12).transpose()) - 4*Conj(Lambdax)*TLambdax*(TLambda12.conjugate()*
      (Lambda12).transpose()) - 2*Conj(Sigmax)*TSigmax*(TLambda12.conjugate()*(
      Lambda12).transpose()) - 6*traceKappaAdjKappa*(TLambda12.conjugate()*(
      TLambda12).transpose()) - 4*traceLambda12AdjLambda12*(TLambda12.conjugate
      ()*(TLambda12).transpose()) - 4*AbsSqr(Lambdax)*(TLambda12.conjugate()*(
      TLambda12).transpose()) - 2*AbsSqr(Sigmax)*(TLambda12.conjugate()*(
      TLambda12).transpose()) + 0.5*Sqr(g1p)*(TLambda12.conjugate()*(TLambda12)
      .transpose()) + 0.1*Sqr(g1p)*Sqr(QS)*(TLambda12.conjugate()*(TLambda12)
      .transpose()) - tracefdAdjfd*(mH2I2*fd.adjoint()*fd) - 3*traceYdAdjYd*(
      mH2I2*fd.adjoint()*fd) - traceYeAdjYe*(mH2I2*fd.adjoint()*fd) - AbsSqr(
      Lambdax)*(mH2I2*fd.adjoint()*fd) + 1.5*Sqr(g1p)*(mH2I2*fd.adjoint()*fd) -
      3*traceKappaAdjKappa*(mH2I2*Lambda12.conjugate()*(Lambda12).transpose())
      - 2*traceLambda12AdjLambda12*(mH2I2*Lambda12.conjugate()*(Lambda12)
      .transpose()) - 2*AbsSqr(Lambdax)*(mH2I2*Lambda12.conjugate()*(Lambda12)
      .transpose()) - AbsSqr(Sigmax)*(mH2I2*Lambda12.conjugate()*(Lambda12)
      .transpose()) + 0.25*Sqr(g1p)*(mH2I2*Lambda12.conjugate()*(Lambda12)
      .transpose()) + 0.05*Sqr(g1p)*Sqr(QS)*(mH2I2*Lambda12.conjugate()*(
      Lambda12).transpose()) - tracefdAdjfd*(fd.adjoint()*fd*mH2I2) - 3*
      traceYdAdjYd*(fd.adjoint()*fd*mH2I2) - traceYeAdjYe*(fd.adjoint()*fd*
      mH2I2) - AbsSqr(Lambdax)*(fd.adjoint()*fd*mH2I2) + 1.5*Sqr(g1p)*(
      fd.adjoint()*fd*mH2I2) - 2*tracefdAdjfd*(fd.adjoint()*mSI2.conjugate()*fd
      ) - 6*traceYdAdjYd*(fd.adjoint()*mSI2.conjugate()*fd) - 2*traceYeAdjYe*(
      fd.adjoint()*mSI2.conjugate()*fd) - 2*AbsSqr(Lambdax)*(fd.adjoint()*
      mSI2.conjugate()*fd) + 3*Sqr(g1p)*(fd.adjoint()*mSI2.conjugate()*fd) - 6*
      traceKappaAdjKappa*(Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12)
      .transpose()) - 4*traceLambda12AdjLambda12*(Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) - 4*AbsSqr(Lambdax)*(
      Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12).transpose()) - 2*AbsSqr
      (Sigmax)*(Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12).transpose())
      + 0.5*Sqr(g1p)*(Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12)
      .transpose()) + 0.1*Sqr(g1p)*Sqr(QS)*(Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) - 3*traceKappaAdjKappa*(
      Lambda12.conjugate()*(Lambda12).transpose()*mH2I2) - 2*
      traceLambda12AdjLambda12*(Lambda12.conjugate()*(Lambda12).transpose()*
      mH2I2) - 2*AbsSqr(Lambdax)*(Lambda12.conjugate()*(Lambda12).transpose()*
      mH2I2) - AbsSqr(Sigmax)*(Lambda12.conjugate()*(Lambda12).transpose()*
      mH2I2) + 0.25*Sqr(g1p)*(Lambda12.conjugate()*(Lambda12).transpose()*mH2I2
      ) + 0.05*Sqr(g1p)*Sqr(QS)*(Lambda12.conjugate()*(Lambda12).transpose()*
      mH2I2) - 8*mHd2*(fd.adjoint()*fd*fd.adjoint()*fd) - 4*(fd.adjoint()*fd*(
      Tfd).adjoint()*Tfd) - 4*mHd2*(fd.adjoint()*fu*fu.adjoint()*fd) - 4*mHu2*(
      fd.adjoint()*fu*fu.adjoint()*fd) - 4*(fd.adjoint()*fu*(Tfu).adjoint()*Tfd
      ) - 4*(fd.adjoint()*Tfd*(Tfd).adjoint()*fd) - 4*(fd.adjoint()*Tfu*(Tfu)
      .adjoint()*fd) - 4*((Tfd).adjoint()*fd*fd.adjoint()*Tfd) - 4*((Tfd)
      .adjoint()*fu*fu.adjoint()*Tfd) - 4*((Tfd).adjoint()*Tfd*fd.adjoint()*fd)
      - 4*((Tfd).adjoint()*Tfu*fu.adjoint()*fd) - 2*mHu2*(Lambda12.conjugate()
      *fu.transpose()*fu.conjugate()*(Lambda12).transpose()) - 2*ms2*(
      Lambda12.conjugate()*fu.transpose()*fu.conjugate()*(Lambda12).transpose()
      ) - 2*(Lambda12.conjugate()*fu.transpose()*Tfu.conjugate()*(TLambda12)
      .transpose()) - 2*mHp2*(Lambda12.conjugate()*hE.transpose()*hE.conjugate(
      )*(Lambda12).transpose()) - 2*ms2*(Lambda12.conjugate()*hE.transpose()*
      hE.conjugate()*(Lambda12).transpose()) - 2*(Lambda12.conjugate()*
      hE.transpose()*ThE.conjugate()*(TLambda12).transpose()) - 4*ms2*(
      Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(
      Lambda12).transpose()) - 2*(Lambda12.conjugate()*(Lambda12).transpose()*
      TLambda12.conjugate()*(TLambda12).transpose()) - 2*(Lambda12.conjugate()*
      (Tfu).transpose()*Tfu.conjugate()*(Lambda12).transpose()) - 2*(
      Lambda12.conjugate()*(ThE).transpose()*ThE.conjugate()*(Lambda12)
      .transpose()) - 2*(Lambda12.conjugate()*(TLambda12).transpose()*
      TLambda12.conjugate()*(Lambda12).transpose()) - 2*(TLambda12.conjugate()*
      fu.transpose()*fu.conjugate()*(TLambda12).transpose()) - 2*(
      TLambda12.conjugate()*hE.transpose()*hE.conjugate()*(TLambda12).transpose
      ()) - 2*(TLambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate(
      )*(TLambda12).transpose()) - 2*(TLambda12.conjugate()*(Tfu).transpose()*
      fu.conjugate()*(Lambda12).transpose()) - 2*(TLambda12.conjugate()*(ThE)
      .transpose()*hE.conjugate()*(Lambda12).transpose()) - 2*(
      TLambda12.conjugate()*(TLambda12).transpose()*Lambda12.conjugate()*(
      Lambda12).transpose()) - 2*(mH2I2*fd.adjoint()*fd*fd.adjoint()*fd) - 2*(
      mH2I2*fd.adjoint()*fu*fu.adjoint()*fd) - mH2I2*Lambda12.conjugate()*
      fu.transpose()*fu.conjugate()*(Lambda12).transpose() - mH2I2*
      Lambda12.conjugate()*hE.transpose()*hE.conjugate()*(Lambda12).transpose()
      - mH2I2*Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()
      *(Lambda12).transpose() - 4*(fd.adjoint()*fd*mH2I2*fd.adjoint()*fd) - 2*(
      fd.adjoint()*fd*fd.adjoint()*fd*mH2I2) - 4*(fd.adjoint()*fd*fd.adjoint()*
      mSI2.conjugate()*fd) - 4*(fd.adjoint()*fu*mH1I2*fu.adjoint()*fd) - 2*(
      fd.adjoint()*fu*fu.adjoint()*fd*mH2I2) - 4*(fd.adjoint()*fu*fu.adjoint()*
      mSI2.conjugate()*fd) - 4*(fd.adjoint()*mSI2.conjugate()*fd*fd.adjoint()*
      fd) - 4*(fd.adjoint()*mSI2.conjugate()*fu*fu.adjoint()*fd) - 2*(
      Lambda12.conjugate()*mH1I2.conjugate()*fu.transpose()*fu.conjugate()*(
      Lambda12).transpose()) - 2*(Lambda12.conjugate()*mH1I2.conjugate()*
      hE.transpose()*hE.conjugate()*(Lambda12).transpose()) - 2*(
      Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12).transpose()*
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*(Lambda12.conjugate()*
      fu.transpose()*mSI2*fu.conjugate()*(Lambda12).transpose()) - 2*(
      Lambda12.conjugate()*fu.transpose()*fu.conjugate()*mH1I2.conjugate()*(
      Lambda12).transpose()) - Lambda12.conjugate()*fu.transpose()*fu.conjugate
      ()*(Lambda12).transpose()*mH2I2 - 2*(Lambda12.conjugate()*hE.transpose()*
      hE.conjugate()*mH1I2.conjugate()*(Lambda12).transpose()) -
      Lambda12.conjugate()*hE.transpose()*hE.conjugate()*(Lambda12).transpose()
      *mH2I2 - 2*(Lambda12.conjugate()*hE.transpose()*me2.conjugate()*
      hE.conjugate()*(Lambda12).transpose()) - 2*(Lambda12.conjugate()*(
      Lambda12).transpose()*mH2I2*Lambda12.conjugate()*(Lambda12).transpose())
      - 2*(Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) - Lambda12.conjugate()*(
      Lambda12).transpose()*Lambda12.conjugate()*(Lambda12).transpose()*mH2I2 +
      6*Power(g2,4)*Tr22*UNITMATRIX(2) - 0.9797958971132712*g1*g1p*Tr2U114*
      UNITMATRIX(2) - 0.9797958971132712*g1*g1p*Tr2U141*UNITMATRIX(2) +
      3.0983866769659336*g1*Tr31*UNITMATRIX(2) - 2.5298221281347035*g1p*Tr34*
      UNITMATRIX(2) + 87*Power(g2,4)*AbsSqr(MassWB)*UNITMATRIX(2) + 1.2*Tr2U111
      *Sqr(g1)*UNITMATRIX(2) + 0.8*Tr2U144*Sqr(g1p)*UNITMATRIX(2) + 3.6*AbsSqr(
      MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*
      Sqr(g2)*UNITMATRIX(2) + 2.4*AbsSqr(MassWB)*Sqr(g1p)*Sqr(g2)*UNITMATRIX(2)
      + 1.2*MassBp*Conj(MassWB)*Sqr(g1p)*Sqr(g2)*UNITMATRIX(2) + 0.36*Conj(
      MassB)*Sqr(g1)*(99*MassB*Sqr(g1) + 2*(2*MassB + MassBp)*Sqr(g1p) + 5*(2*
      MassB + MassWB)*Sqr(g2))*UNITMATRIX(2) + 0.02*Conj(MassBp)*Sqr(g1p)*(5*(
      60*MassBp*(fd.adjoint()*fd) - 30*(fd.adjoint()*Tfd) + (5 + Sqr(QS))*(2*
      MassBp*(Lambda12.conjugate()*(Lambda12).transpose()) - Lambda12.conjugate
      ()*(TLambda12).transpose())) + 6*(6*(MassB + 2*MassBp)*Sqr(g1) + 10*(2*
      MassBp + MassWB)*Sqr(g2) + MassBp*Sqr(g1p)*(192 + Sqr(QS)))*UNITMATRIX(2)
      ))).real();


   return beta_mH2I2;
}

} // namespace flexiblesusy
