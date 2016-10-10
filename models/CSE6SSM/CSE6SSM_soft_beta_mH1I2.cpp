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

// File generated at Wed 3 Jun 2015 23:43:33

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mH1I2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_mH1I2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = (oneOver16PiSqr*(2*mHu2*(fu.adjoint()*fu) + 2*mHp2*(
      hE.adjoint()*hE) + 2*ms2*((Lambda12).adjoint()*Lambda12) + 2*((Tfu)
      .adjoint()*Tfu) + 2*((ThE).adjoint()*ThE) + 2*((TLambda12).adjoint()*
      TLambda12) + mH1I2*fu.adjoint()*fu + mH1I2*hE.adjoint()*hE + mH1I2*(
      Lambda12).adjoint()*Lambda12 + fu.adjoint()*fu*mH1I2 + 2*(fu.adjoint()*
      mSI2.conjugate()*fu) + hE.adjoint()*hE*mH1I2 + 2*(hE.adjoint()*me2*hE) +
      2*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) + (Lambda12).adjoint(
      )*Lambda12*mH1I2 - 0.7745966692414834*g1*Tr11*UNITMATRIX(2) -
      0.9486832980505138*g1p*Tr14*UNITMATRIX(2) - 1.2*AbsSqr(MassB)*Sqr(g1)*
      UNITMATRIX(2) - 1.8*AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(2) - 6*AbsSqr(
      MassWB)*Sqr(g2)*UNITMATRIX(2))).real();


   return beta_mH1I2;
}

/**
 * Calculates the two-loop beta function of mH1I2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CSE6SSM_soft_parameters::calc_beta_mH1I2_two_loop(const Soft_traces& soft_traces) const
{
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
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDAdjgDconjmq2 = TRACE_STRUCT.tracegDAdjgDconjmq2;
   const double tracegDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgD;
   const double tracehEmH1I2AdjhE = TRACE_STRUCT.tracehEmH1I2AdjhE;
   const double tracehEAdjhEme2 = TRACE_STRUCT.tracehEAdjhEme2;
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
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
   const double traceAdjfuTfu = TRACE_STRUCT.traceAdjfuTfu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
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


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = (twoLoop*(-2*traceconjTfuTpTfu*(fu.adjoint()*fu) - 6*
      traceconjTYuTpTYu*(fu.adjoint()*fu) - 4*mHu2*tracefuAdjfu*(fu.adjoint()*
      fu) - 2*tracefuAdjfuconjmSI2*(fu.adjoint()*fu) - 2*tracefumH1I2Adjfu*(
      fu.adjoint()*fu) - 6*tracemq2AdjYuYu*(fu.adjoint()*fu) - 6*
      tracemu2YuAdjYu*(fu.adjoint()*fu) - 12*mHu2*traceYuAdjYu*(fu.adjoint()*fu
      ) - 2*mHd2*AbsSqr(Lambdax)*(fu.adjoint()*fu) - 4*mHu2*AbsSqr(Lambdax)*(
      fu.adjoint()*fu) - 2*ms2*AbsSqr(Lambdax)*(fu.adjoint()*fu) - 2*AbsSqr(
      TLambdax)*(fu.adjoint()*fu) + 2*mHu2*Sqr(g1p)*(fu.adjoint()*fu) - 2*
      traceconjTfuTpfu*(fu.adjoint()*Tfu) - 6*traceconjTYuTpYu*(fu.adjoint()*
      Tfu) - 2*Conj(TLambdax)*Lambdax*(fu.adjoint()*Tfu) - 6*traceconjTgDTpTgD*
      (hE.adjoint()*hE) - 2*traceconjThETpThE*(hE.adjoint()*hE) - 12*mHp2*
      tracegDAdjgD*(hE.adjoint()*hE) - 6*tracegDAdjgDconjmq2*(hE.adjoint()*hE)
      - 6*tracegDconjmDxbar2AdjgD*(hE.adjoint()*hE) - 4*mHp2*tracehEAdjhE*(
      hE.adjoint()*hE) - 2*tracehEAdjhEme2*(hE.adjoint()*hE) - 2*
      tracehEmH1I2AdjhE*(hE.adjoint()*hE) - 4*mHp2*AbsSqr(SigmaL)*(hE.adjoint()
      *hE) - 2*mHpbar2*AbsSqr(SigmaL)*(hE.adjoint()*hE) - 2*mphi2*AbsSqr(SigmaL
      )*(hE.adjoint()*hE) - 2*AbsSqr(TSigmaL)*(hE.adjoint()*hE) + 2.4*mHp2*Sqr(
      g1)*(hE.adjoint()*hE) - 0.4*mHp2*Sqr(g1p)*(hE.adjoint()*hE) - 6*
      traceconjTgDTpgD*(hE.adjoint()*ThE) - 2*traceconjThETphE*(hE.adjoint()*
      ThE) - 2*Conj(TSigmaL)*SigmaL*(hE.adjoint()*ThE) - 6*
      traceconjTKappaTpTKappa*((Lambda12).adjoint()*Lambda12) - 4*
      traceconjTLambda12TpTLambda12*((Lambda12).adjoint()*Lambda12) - 12*ms2*
      traceKappaAdjKappa*((Lambda12).adjoint()*Lambda12) - 6*
      traceKappaAdjKappaconjmDx2*((Lambda12).adjoint()*Lambda12) - 6*
      traceKappaconjmDxbar2AdjKappa*((Lambda12).adjoint()*Lambda12) - 8*ms2*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*Lambda12) - 4*
      traceLambda12AdjLambda12conjmH2I2*((Lambda12).adjoint()*Lambda12) - 4*
      tracemH1I2AdjLambda12Lambda12*((Lambda12).adjoint()*Lambda12) - 4*mHd2*
      AbsSqr(Lambdax)*((Lambda12).adjoint()*Lambda12) - 4*mHu2*AbsSqr(Lambdax)*
      ((Lambda12).adjoint()*Lambda12) - 8*ms2*AbsSqr(Lambdax)*((Lambda12)
      .adjoint()*Lambda12) - 2*mphi2*AbsSqr(Sigmax)*((Lambda12).adjoint()*
      Lambda12) - 4*ms2*AbsSqr(Sigmax)*((Lambda12).adjoint()*Lambda12) - 2*
      msbar2*AbsSqr(Sigmax)*((Lambda12).adjoint()*Lambda12) - 4*AbsSqr(TLambdax
      )*((Lambda12).adjoint()*Lambda12) - 2*AbsSqr(TSigmax)*((Lambda12).adjoint
      ()*Lambda12) - 0.5*ms2*Sqr(g1p)*((Lambda12).adjoint()*Lambda12) + 0.1*ms2
      *Sqr(g1p)*Sqr(QS)*((Lambda12).adjoint()*Lambda12) - 6*
      traceconjTKappaTpKappa*((Lambda12).adjoint()*TLambda12) - 4*
      traceconjTLambda12TpLambda12*((Lambda12).adjoint()*TLambda12) - 4*Conj(
      TLambdax)*Lambdax*((Lambda12).adjoint()*TLambda12) - 2*Conj(TSigmax)*
      Sigmax*((Lambda12).adjoint()*TLambda12) - 2*traceAdjfuTfu*((Tfu).adjoint(
      )*fu) - 6*traceAdjYuTYu*((Tfu).adjoint()*fu) - 2*MassBp*Sqr(g1p)*((Tfu)
      .adjoint()*fu) - 2*Conj(Lambdax)*TLambdax*((Tfu).adjoint()*fu) - 2*
      tracefuAdjfu*((Tfu).adjoint()*Tfu) - 6*traceYuAdjYu*((Tfu).adjoint()*Tfu)
      - 2*AbsSqr(Lambdax)*((Tfu).adjoint()*Tfu) + 2*Sqr(g1p)*((Tfu).adjoint()*
      Tfu) - 6*traceAdjgDTgD*((ThE).adjoint()*hE) - 2*traceAdjhEThE*((ThE)
      .adjoint()*hE) - 2.4*MassB*Sqr(g1)*((ThE).adjoint()*hE) + 0.4*MassBp*Sqr(
      g1p)*((ThE).adjoint()*hE) - 2*Conj(SigmaL)*TSigmaL*((ThE).adjoint()*hE) -
      6*tracegDAdjgD*((ThE).adjoint()*ThE) - 2*tracehEAdjhE*((ThE).adjoint()*
      ThE) - 2*AbsSqr(SigmaL)*((ThE).adjoint()*ThE) + 2.4*Sqr(g1)*((ThE)
      .adjoint()*ThE) - 0.4*Sqr(g1p)*((ThE).adjoint()*ThE) - 6*
      traceAdjKappaTKappa*((TLambda12).adjoint()*Lambda12) - 4*
      traceAdjLambda12TLambda12*((TLambda12).adjoint()*Lambda12) + 0.5*MassBp*
      Sqr(g1p)*((TLambda12).adjoint()*Lambda12) - 0.1*MassBp*Sqr(g1p)*Sqr(QS)*(
      (TLambda12).adjoint()*Lambda12) - 4*Conj(Lambdax)*TLambdax*((TLambda12)
      .adjoint()*Lambda12) - 2*Conj(Sigmax)*TSigmax*((TLambda12).adjoint()*
      Lambda12) - 6*traceKappaAdjKappa*((TLambda12).adjoint()*TLambda12) - 4*
      traceLambda12AdjLambda12*((TLambda12).adjoint()*TLambda12) - 4*AbsSqr(
      Lambdax)*((TLambda12).adjoint()*TLambda12) - 2*AbsSqr(Sigmax)*((TLambda12
      ).adjoint()*TLambda12) - 0.5*Sqr(g1p)*((TLambda12).adjoint()*TLambda12) +
      0.1*Sqr(g1p)*Sqr(QS)*((TLambda12).adjoint()*TLambda12) - tracefuAdjfu*(
      mH1I2*fu.adjoint()*fu) - 3*traceYuAdjYu*(mH1I2*fu.adjoint()*fu) - AbsSqr(
      Lambdax)*(mH1I2*fu.adjoint()*fu) + Sqr(g1p)*(mH1I2*fu.adjoint()*fu) - 3*
      tracegDAdjgD*(mH1I2*hE.adjoint()*hE) - tracehEAdjhE*(mH1I2*hE.adjoint()*
      hE) - AbsSqr(SigmaL)*(mH1I2*hE.adjoint()*hE) + 1.2*Sqr(g1)*(mH1I2*
      hE.adjoint()*hE) - 0.2*Sqr(g1p)*(mH1I2*hE.adjoint()*hE) - 3*
      traceKappaAdjKappa*(mH1I2*(Lambda12).adjoint()*Lambda12) - 2*
      traceLambda12AdjLambda12*(mH1I2*(Lambda12).adjoint()*Lambda12) - 2*AbsSqr
      (Lambdax)*(mH1I2*(Lambda12).adjoint()*Lambda12) - AbsSqr(Sigmax)*(mH1I2*(
      Lambda12).adjoint()*Lambda12) - 0.25*Sqr(g1p)*(mH1I2*(Lambda12).adjoint()
      *Lambda12) + 0.05*Sqr(g1p)*Sqr(QS)*(mH1I2*(Lambda12).adjoint()*Lambda12)
      - tracefuAdjfu*(fu.adjoint()*fu*mH1I2) - 3*traceYuAdjYu*(fu.adjoint()*fu*
      mH1I2) - AbsSqr(Lambdax)*(fu.adjoint()*fu*mH1I2) + Sqr(g1p)*(fu.adjoint()
      *fu*mH1I2) - 2*tracefuAdjfu*(fu.adjoint()*mSI2.conjugate()*fu) - 6*
      traceYuAdjYu*(fu.adjoint()*mSI2.conjugate()*fu) - 2*AbsSqr(Lambdax)*(
      fu.adjoint()*mSI2.conjugate()*fu) + 2*Sqr(g1p)*(fu.adjoint()*
      mSI2.conjugate()*fu) - 3*tracegDAdjgD*(hE.adjoint()*hE*mH1I2) -
      tracehEAdjhE*(hE.adjoint()*hE*mH1I2) - AbsSqr(SigmaL)*(hE.adjoint()*hE*
      mH1I2) + 1.2*Sqr(g1)*(hE.adjoint()*hE*mH1I2) - 0.2*Sqr(g1p)*(hE.adjoint()
      *hE*mH1I2) - 6*tracegDAdjgD*(hE.adjoint()*me2*hE) - 2*tracehEAdjhE*(
      hE.adjoint()*me2*hE) - 2*AbsSqr(SigmaL)*(hE.adjoint()*me2*hE) + 2.4*Sqr(
      g1)*(hE.adjoint()*me2*hE) - 0.4*Sqr(g1p)*(hE.adjoint()*me2*hE) - 6*
      traceKappaAdjKappa*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - 4*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12
      ) - 4*AbsSqr(Lambdax)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) -
      2*AbsSqr(Sigmax)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - 0.5
      *Sqr(g1p)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) + 0.1*Sqr(g1p
      )*Sqr(QS)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - 3*
      traceKappaAdjKappa*((Lambda12).adjoint()*Lambda12*mH1I2) - 2*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*Lambda12*mH1I2) - 2*AbsSqr
      (Lambdax)*((Lambda12).adjoint()*Lambda12*mH1I2) - AbsSqr(Sigmax)*((
      Lambda12).adjoint()*Lambda12*mH1I2) - 0.25*Sqr(g1p)*((Lambda12).adjoint()
      *Lambda12*mH1I2) + 0.05*Sqr(g1p)*Sqr(QS)*((Lambda12).adjoint()*Lambda12*
      mH1I2) - 4*mHd2*(fu.adjoint()*fd*fd.adjoint()*fu) - 4*mHu2*(fu.adjoint()*
      fd*fd.adjoint()*fu) - 4*(fu.adjoint()*fd*(Tfd).adjoint()*Tfu) - 8*mHu2*(
      fu.adjoint()*fu*fu.adjoint()*fu) - 4*(fu.adjoint()*fu*(Tfu).adjoint()*Tfu
      ) - 4*(fu.adjoint()*Tfd*(Tfd).adjoint()*fu) - 4*(fu.adjoint()*Tfu*(Tfu)
      .adjoint()*fu) - 8*mHp2*(hE.adjoint()*hE*hE.adjoint()*hE) - 4*(hE.adjoint
      ()*hE*(ThE).adjoint()*ThE) - 4*mHd2*(hE.adjoint()*Ye*Ye.adjoint()*hE) - 4
      *mHp2*(hE.adjoint()*Ye*Ye.adjoint()*hE) - 4*(hE.adjoint()*Ye*(TYe)
      .adjoint()*ThE) - 4*(hE.adjoint()*ThE*(ThE).adjoint()*hE) - 4*(hE.adjoint
      ()*TYe*(TYe).adjoint()*hE) - 4*ms2*((Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*Lambda12*(
      TLambda12).adjoint()*TLambda12) - 2*((Lambda12).adjoint()*TLambda12*(
      TLambda12).adjoint()*Lambda12) - 2*mHd2*((Lambda12).adjoint()*
      fd.transpose()*fd.conjugate()*Lambda12) - 2*ms2*((Lambda12).adjoint()*
      fd.transpose()*fd.conjugate()*Lambda12) - 2*((Lambda12).adjoint()*
      fd.transpose()*Tfd.conjugate()*TLambda12) - 2*((Lambda12).adjoint()*(Tfd)
      .transpose()*Tfd.conjugate()*Lambda12) - 4*((Tfu).adjoint()*fd*fd.adjoint
      ()*Tfu) - 4*((Tfu).adjoint()*fu*fu.adjoint()*Tfu) - 4*((Tfu).adjoint()*
      Tfd*fd.adjoint()*fu) - 4*((Tfu).adjoint()*Tfu*fu.adjoint()*fu) - 4*((ThE)
      .adjoint()*hE*hE.adjoint()*ThE) - 4*((ThE).adjoint()*Ye*Ye.adjoint()*ThE)
      - 4*((ThE).adjoint()*ThE*hE.adjoint()*hE) - 4*((ThE).adjoint()*TYe*
      Ye.adjoint()*hE) - 2*((TLambda12).adjoint()*Lambda12*(Lambda12).adjoint()
      *TLambda12) - 2*((TLambda12).adjoint()*TLambda12*(Lambda12).adjoint()*
      Lambda12) - 2*((TLambda12).adjoint()*fd.transpose()*fd.conjugate()*
      TLambda12) - 2*((TLambda12).adjoint()*(Tfd).transpose()*fd.conjugate()*
      Lambda12) - 2*(mH1I2*fu.adjoint()*fd*fd.adjoint()*fu) - 2*(mH1I2*
      fu.adjoint()*fu*fu.adjoint()*fu) - 2*(mH1I2*hE.adjoint()*hE*hE.adjoint()*
      hE) - 2*(mH1I2*hE.adjoint()*Ye*Ye.adjoint()*hE) - mH1I2*(Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12 - mH1I2*(Lambda12)
      .adjoint()*fd.transpose()*fd.conjugate()*Lambda12 - 4*(fu.adjoint()*fd*
      mH2I2*fd.adjoint()*fu) - 2*(fu.adjoint()*fd*fd.adjoint()*fu*mH1I2) - 4*(
      fu.adjoint()*fd*fd.adjoint()*mSI2.conjugate()*fu) - 4*(fu.adjoint()*fu*
      mH1I2*fu.adjoint()*fu) - 2*(fu.adjoint()*fu*fu.adjoint()*fu*mH1I2) - 4*(
      fu.adjoint()*fu*fu.adjoint()*mSI2.conjugate()*fu) - 4*(fu.adjoint()*
      mSI2.conjugate()*fd*fd.adjoint()*fu) - 4*(fu.adjoint()*mSI2.conjugate()*
      fu*fu.adjoint()*fu) - 4*(hE.adjoint()*hE*mH1I2*hE.adjoint()*hE) - 2*(
      hE.adjoint()*hE*hE.adjoint()*hE*mH1I2) - 4*(hE.adjoint()*hE*hE.adjoint()*
      me2*hE) - 4*(hE.adjoint()*me2*hE*hE.adjoint()*hE) - 4*(hE.adjoint()*me2*
      Ye*Ye.adjoint()*hE) - 4*(hE.adjoint()*Ye*ml2*Ye.adjoint()*hE) - 2*(
      hE.adjoint()*Ye*Ye.adjoint()*hE*mH1I2) - 4*(hE.adjoint()*Ye*Ye.adjoint()*
      me2*hE) - 2*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12*(Lambda12)
      .adjoint()*Lambda12) - 2*((Lambda12).adjoint()*mH2I2.conjugate()*
      fd.transpose()*fd.conjugate()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*mH1I2*(Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*(Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - (Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12*mH1I2 - 2*((Lambda12)
      .adjoint()*fd.transpose()*mSI2*fd.conjugate()*Lambda12) - 2*((Lambda12)
      .adjoint()*fd.transpose()*fd.conjugate()*mH2I2.conjugate()*Lambda12) - (
      Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12*mH1I2 + 6*
      Power(g2,4)*Tr22*UNITMATRIX(2) + 1.4696938456699067*g1*g1p*Tr2U114*
      UNITMATRIX(2) + 1.4696938456699067*g1*g1p*Tr2U141*UNITMATRIX(2) -
      3.0983866769659336*g1*Tr31*UNITMATRIX(2) - 3.794733192202055*g1p*Tr34*
      UNITMATRIX(2) + 87*Power(g2,4)*AbsSqr(MassWB)*UNITMATRIX(2) + 1.2*Tr2U111
      *Sqr(g1)*UNITMATRIX(2) + 1.8*Tr2U144*Sqr(g1p)*UNITMATRIX(2) + 3.6*AbsSqr(
      MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*
      Sqr(g2)*UNITMATRIX(2) + 5.4*AbsSqr(MassWB)*Sqr(g1p)*Sqr(g2)*UNITMATRIX(2)
      + 2.7*MassBp*Conj(MassWB)*Sqr(g1p)*Sqr(g2)*UNITMATRIX(2) + 0.06*Conj(
      MassB)*Sqr(g1)*(80*MassB*(hE.adjoint()*hE) - 40*(hE.adjoint()*ThE) + (594
      *MassB*Sqr(g1) - 3*(2*MassB + MassBp)*Sqr(g1p) + 30*(2*MassB + MassWB)*
      Sqr(g2))*UNITMATRIX(2)) + 0.01*Conj(MassBp)*Sqr(g1p)*(10*(40*MassBp*(
      fu.adjoint()*fu) - 20*(fu.adjoint()*Tfu) - 8*MassBp*(hE.adjoint()*hE) + 4
      *(hE.adjoint()*ThE) - 10*MassBp*((Lambda12).adjoint()*Lambda12) + 2*
      MassBp*Sqr(QS)*((Lambda12).adjoint()*Lambda12) + 5*((Lambda12).adjoint()*
      TLambda12) - Sqr(QS)*((Lambda12).adjoint()*TLambda12)) - 9*(2*(MassB + 2*
      MassBp)*Sqr(g1) - 3*(10*(2*MassBp + MassWB)*Sqr(g2) + MassBp*Sqr(g1p)*(
      197 + Sqr(QS))))*UNITMATRIX(2)))).real();


   return beta_mH1I2;
}

} // namespace flexiblesusy
