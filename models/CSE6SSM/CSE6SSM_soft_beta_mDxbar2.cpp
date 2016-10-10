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

// File generated at Wed 3 Jun 2015 23:43:39

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mDxbar2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mDxbar2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = (oneOver16PiSqr*(4*mHp2*(gD.transpose()*gD.conjugate())
      + 2*ms2*((Kappa).transpose()*Kappa.conjugate()) + 4*((TgD).transpose()*
      TgD.conjugate()) + 2*((TKappa).transpose()*TKappa.conjugate()) + 2*(
      mDxbar2*gD.transpose()*gD.conjugate()) + mDxbar2*(Kappa).transpose()*
      Kappa.conjugate() + 4*(gD.transpose()*mq2*gD.conjugate()) + 2*(
      gD.transpose()*gD.conjugate()*mDxbar2) + 2*((Kappa).transpose()*mDx2*
      Kappa.conjugate()) + (Kappa).transpose()*Kappa.conjugate()*mDxbar2 +
      0.5163977794943222*g1*Tr11*UNITMATRIX(3) - 0.9486832980505138*g1p*Tr14*
      UNITMATRIX(3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) -
      1.8*AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(
      MassG)*Sqr(g3)*UNITMATRIX(3))).real();


   return beta_mDxbar2;
}

/**
 * Calculates the two-loop beta function of mDxbar2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_mDxbar2_two_loop(const Soft_traces& soft_traces) const
{
   const double tracegDAdjgD = TRACE_STRUCT.tracegDAdjgD;
   const double tracehEAdjhE = TRACE_STRUCT.tracehEAdjhE;
   const double traceconjTgDTpTgD = TRACE_STRUCT.traceconjTgDTpTgD;
   const double traceconjThETpThE = TRACE_STRUCT.traceconjThETpThE;
   const double tracegDAdjgDconjmq2 = TRACE_STRUCT.tracegDAdjgDconjmq2;
   const double tracegDconjmDxbar2AdjgD =
      TRACE_STRUCT.tracegDconjmDxbar2AdjgD;
   const double tracehEmH1I2AdjhE = TRACE_STRUCT.tracehEmH1I2AdjhE;
   const double tracehEAdjhEme2 = TRACE_STRUCT.tracehEAdjhEme2;
   const double traceAdjgDTgD = TRACE_STRUCT.traceAdjgDTgD;
   const double traceAdjhEThE = TRACE_STRUCT.traceAdjhEThE;
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
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = (twoLoop*(-12*traceconjTgDTpTgD*(gD.transpose()*
      gD.conjugate()) - 4*traceconjThETpThE*(gD.transpose()*gD.conjugate()) -
      24*mHp2*tracegDAdjgD*(gD.transpose()*gD.conjugate()) - 12*
      tracegDAdjgDconjmq2*(gD.transpose()*gD.conjugate()) - 12*
      tracegDconjmDxbar2AdjgD*(gD.transpose()*gD.conjugate()) - 8*mHp2*
      tracehEAdjhE*(gD.transpose()*gD.conjugate()) - 4*tracehEAdjhEme2*(
      gD.transpose()*gD.conjugate()) - 4*tracehEmH1I2AdjhE*(gD.transpose()*
      gD.conjugate()) - 8*mHp2*AbsSqr(SigmaL)*(gD.transpose()*gD.conjugate()) -
      4*mHpbar2*AbsSqr(SigmaL)*(gD.transpose()*gD.conjugate()) - 4*mphi2*
      AbsSqr(SigmaL)*(gD.transpose()*gD.conjugate()) - 4*AbsSqr(TSigmaL)*(
      gD.transpose()*gD.conjugate()) + 0.8*mHp2*Sqr(g1)*(gD.transpose()*
      gD.conjugate()) - 0.8*mHp2*Sqr(g1p)*(gD.transpose()*gD.conjugate()) + 12*
      mHp2*Sqr(g2)*(gD.transpose()*gD.conjugate()) + 24*AbsSqr(MassWB)*Sqr(g2)*
      (gD.transpose()*gD.conjugate()) - 12*traceAdjgDTgD*(gD.transpose()*
      TgD.conjugate()) - 4*traceAdjhEThE*(gD.transpose()*TgD.conjugate()) - 0.8
      *MassB*Sqr(g1)*(gD.transpose()*TgD.conjugate()) + 0.8*MassBp*Sqr(g1p)*(
      gD.transpose()*TgD.conjugate()) - 12*MassWB*Sqr(g2)*(gD.transpose()*
      TgD.conjugate()) - 4*Conj(SigmaL)*TSigmaL*(gD.transpose()*TgD.conjugate()
      ) - 6*traceconjTKappaTpTKappa*((Kappa).transpose()*Kappa.conjugate()) - 4
      *traceconjTLambda12TpTLambda12*((Kappa).transpose()*Kappa.conjugate()) -
      12*ms2*traceKappaAdjKappa*((Kappa).transpose()*Kappa.conjugate()) - 6*
      traceKappaAdjKappaconjmDx2*((Kappa).transpose()*Kappa.conjugate()) - 6*
      traceKappaconjmDxbar2AdjKappa*((Kappa).transpose()*Kappa.conjugate()) - 8
      *ms2*traceLambda12AdjLambda12*((Kappa).transpose()*Kappa.conjugate()) - 4
      *traceLambda12AdjLambda12conjmH2I2*((Kappa).transpose()*Kappa.conjugate()
      ) - 4*tracemH1I2AdjLambda12Lambda12*((Kappa).transpose()*Kappa.conjugate(
      )) - 4*mHd2*AbsSqr(Lambdax)*((Kappa).transpose()*Kappa.conjugate()) - 4*
      mHu2*AbsSqr(Lambdax)*((Kappa).transpose()*Kappa.conjugate()) - 8*ms2*
      AbsSqr(Lambdax)*((Kappa).transpose()*Kappa.conjugate()) - 2*mphi2*AbsSqr(
      Sigmax)*((Kappa).transpose()*Kappa.conjugate()) - 4*ms2*AbsSqr(Sigmax)*((
      Kappa).transpose()*Kappa.conjugate()) - 2*msbar2*AbsSqr(Sigmax)*((Kappa)
      .transpose()*Kappa.conjugate()) - 4*AbsSqr(TLambdax)*((Kappa).transpose()
      *Kappa.conjugate()) - 2*AbsSqr(TSigmax)*((Kappa).transpose()*
      Kappa.conjugate()) - 0.5*ms2*Sqr(g1p)*((Kappa).transpose()*
      Kappa.conjugate()) + 0.1*ms2*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceAdjKappaTKappa*((Kappa).transpose()*
      TKappa.conjugate()) - 4*traceAdjLambda12TLambda12*((Kappa).transpose()*
      TKappa.conjugate()) + 0.5*MassBp*Sqr(g1p)*((Kappa).transpose()*
      TKappa.conjugate()) - 0.1*MassBp*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*
      TKappa.conjugate()) - 4*Conj(Lambdax)*TLambdax*((Kappa).transpose()*
      TKappa.conjugate()) - 2*Conj(Sigmax)*TSigmax*((Kappa).transpose()*
      TKappa.conjugate()) - 12*traceconjTgDTpgD*((TgD).transpose()*gD.conjugate
      ()) - 4*traceconjThETphE*((TgD).transpose()*gD.conjugate()) - 4*Conj(
      TSigmaL)*SigmaL*((TgD).transpose()*gD.conjugate()) - 12*Conj(MassWB)*Sqr(
      g2)*((TgD).transpose()*gD.conjugate()) - 12*tracegDAdjgD*((TgD).transpose
      ()*TgD.conjugate()) - 4*tracehEAdjhE*((TgD).transpose()*TgD.conjugate())
      - 4*AbsSqr(SigmaL)*((TgD).transpose()*TgD.conjugate()) + 0.8*Sqr(g1)*((
      TgD).transpose()*TgD.conjugate()) - 0.8*Sqr(g1p)*((TgD).transpose()*
      TgD.conjugate()) + 12*Sqr(g2)*((TgD).transpose()*TgD.conjugate()) - 6*
      traceconjTKappaTpKappa*((TKappa).transpose()*Kappa.conjugate()) - 4*
      traceconjTLambda12TpLambda12*((TKappa).transpose()*Kappa.conjugate()) - 4
      *Conj(TLambdax)*Lambdax*((TKappa).transpose()*Kappa.conjugate()) - 2*Conj
      (TSigmax)*Sigmax*((TKappa).transpose()*Kappa.conjugate()) - 6*
      traceKappaAdjKappa*((TKappa).transpose()*TKappa.conjugate()) - 4*
      traceLambda12AdjLambda12*((TKappa).transpose()*TKappa.conjugate()) - 4*
      AbsSqr(Lambdax)*((TKappa).transpose()*TKappa.conjugate()) - 2*AbsSqr(
      Sigmax)*((TKappa).transpose()*TKappa.conjugate()) - 0.5*Sqr(g1p)*((TKappa
      ).transpose()*TKappa.conjugate()) + 0.1*Sqr(g1p)*Sqr(QS)*((TKappa)
      .transpose()*TKappa.conjugate()) - 6*tracegDAdjgD*(mDxbar2*gD.transpose()
      *gD.conjugate()) - 2*tracehEAdjhE*(mDxbar2*gD.transpose()*gD.conjugate())
      - 2*AbsSqr(SigmaL)*(mDxbar2*gD.transpose()*gD.conjugate()) + 0.4*Sqr(g1)
      *(mDxbar2*gD.transpose()*gD.conjugate()) - 0.4*Sqr(g1p)*(mDxbar2*
      gD.transpose()*gD.conjugate()) + 6*Sqr(g2)*(mDxbar2*gD.transpose()*
      gD.conjugate()) - 3*traceKappaAdjKappa*(mDxbar2*(Kappa).transpose()*
      Kappa.conjugate()) - 2*traceLambda12AdjLambda12*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 2*AbsSqr(Lambdax)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - AbsSqr(Sigmax)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 0.25*Sqr(g1p)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) + 0.05*Sqr(g1p)*Sqr(QS)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 12*tracegDAdjgD*(gD.transpose()*mq2*
      gD.conjugate()) - 4*tracehEAdjhE*(gD.transpose()*mq2*gD.conjugate()) - 4*
      AbsSqr(SigmaL)*(gD.transpose()*mq2*gD.conjugate()) + 0.8*Sqr(g1)*(
      gD.transpose()*mq2*gD.conjugate()) - 0.8*Sqr(g1p)*(gD.transpose()*mq2*
      gD.conjugate()) + 12*Sqr(g2)*(gD.transpose()*mq2*gD.conjugate()) - 6*
      tracegDAdjgD*(gD.transpose()*gD.conjugate()*mDxbar2) - 2*tracehEAdjhE*(
      gD.transpose()*gD.conjugate()*mDxbar2) - 2*AbsSqr(SigmaL)*(gD.transpose()
      *gD.conjugate()*mDxbar2) + 0.4*Sqr(g1)*(gD.transpose()*gD.conjugate()*
      mDxbar2) - 0.4*Sqr(g1p)*(gD.transpose()*gD.conjugate()*mDxbar2) + 6*Sqr(
      g2)*(gD.transpose()*gD.conjugate()*mDxbar2) - 6*traceKappaAdjKappa*((
      Kappa).transpose()*mDx2*Kappa.conjugate()) - 4*traceLambda12AdjLambda12*(
      (Kappa).transpose()*mDx2*Kappa.conjugate()) - 4*AbsSqr(Lambdax)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 2*AbsSqr(Sigmax)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 0.5*Sqr(g1p)*((Kappa).transpose()*
      mDx2*Kappa.conjugate()) + 0.1*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*mDx2*
      Kappa.conjugate()) - 3*traceKappaAdjKappa*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) - 2*traceLambda12AdjLambda12*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 2*AbsSqr(Lambdax)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - AbsSqr(Sigmax)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 0.25*Sqr(g1p)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) + 0.05*Sqr(g1p)*Sqr(QS)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 4*mHd2*(gD.transpose()*
      Yd.adjoint()*Yd*gD.conjugate()) - 4*mHp2*(gD.transpose()*Yd.adjoint()*Yd*
      gD.conjugate()) - 4*(gD.transpose()*Yd.adjoint()*TYd*TgD.conjugate()) - 4
      *mHp2*(gD.transpose()*Yu.adjoint()*Yu*gD.conjugate()) - 4*mHu2*(
      gD.transpose()*Yu.adjoint()*Yu*gD.conjugate()) - 4*(gD.transpose()*
      Yu.adjoint()*TYu*TgD.conjugate()) - 4*(gD.transpose()*(TYd).adjoint()*TYd
      *gD.conjugate()) - 4*(gD.transpose()*(TYu).adjoint()*TYu*gD.conjugate())
      - 8*mHp2*(gD.transpose()*gD.conjugate()*gD.transpose()*gD.conjugate()) -
      4*(gD.transpose()*gD.conjugate()*(TgD).transpose()*TgD.conjugate()) - 4*(
      gD.transpose()*TgD.conjugate()*(TgD).transpose()*gD.conjugate()) - 4*ms2*
      ((Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()) - 2*((Kappa).transpose()*Kappa.conjugate()*(TKappa)
      .transpose()*TKappa.conjugate()) - 2*((Kappa).transpose()*
      TKappa.conjugate()*(TKappa).transpose()*Kappa.conjugate()) - 4*((TgD)
      .transpose()*Yd.adjoint()*Yd*TgD.conjugate()) - 4*((TgD).transpose()*
      Yu.adjoint()*Yu*TgD.conjugate()) - 4*((TgD).transpose()*(TYd).adjoint()*
      Yd*gD.conjugate()) - 4*((TgD).transpose()*(TYu).adjoint()*Yu*gD.conjugate
      ()) - 4*((TgD).transpose()*gD.conjugate()*gD.transpose()*TgD.conjugate())
      - 4*((TgD).transpose()*TgD.conjugate()*gD.transpose()*gD.conjugate()) -
      2*((TKappa).transpose()*Kappa.conjugate()*(Kappa).transpose()*
      TKappa.conjugate()) - 2*((TKappa).transpose()*TKappa.conjugate()*(Kappa)
      .transpose()*Kappa.conjugate()) - 2*(mDxbar2*gD.transpose()*Yd.adjoint()*
      Yd*gD.conjugate()) - 2*(mDxbar2*gD.transpose()*Yu.adjoint()*Yu*
      gD.conjugate()) - 2*(mDxbar2*gD.transpose()*gD.conjugate()*gD.transpose()
      *gD.conjugate()) - mDxbar2*(Kappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()*Kappa.conjugate() - 4*(gD.transpose()*mq2*Yd.adjoint()*Yd*
      gD.conjugate()) - 4*(gD.transpose()*mq2*Yu.adjoint()*Yu*gD.conjugate()) -
      4*(gD.transpose()*mq2*gD.conjugate()*gD.transpose()*gD.conjugate()) - 4*
      (gD.transpose()*Yd.adjoint()*md2*Yd*gD.conjugate()) - 4*(gD.transpose()*
      Yd.adjoint()*Yd*mq2*gD.conjugate()) - 2*(gD.transpose()*Yd.adjoint()*Yd*
      gD.conjugate()*mDxbar2) - 4*(gD.transpose()*Yu.adjoint()*mu2*Yu*
      gD.conjugate()) - 4*(gD.transpose()*Yu.adjoint()*Yu*mq2*gD.conjugate()) -
      2*(gD.transpose()*Yu.adjoint()*Yu*gD.conjugate()*mDxbar2) - 4*(
      gD.transpose()*gD.conjugate()*mDxbar2*gD.transpose()*gD.conjugate()) - 4*
      (gD.transpose()*gD.conjugate()*gD.transpose()*mq2*gD.conjugate()) - 2*(
      gD.transpose()*gD.conjugate()*gD.transpose()*gD.conjugate()*mDxbar2) - 2*
      ((Kappa).transpose()*mDx2*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()) - 2*((Kappa).transpose()*Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*mDx2*Kappa.conjugate()) - (Kappa)
      .transpose()*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*
      mDxbar2 + 10.666666666666666*Power(g3,4)*Tr23*UNITMATRIX(3) -
      0.9797958971132712*g1*g1p*Tr2U114*UNITMATRIX(3) - 0.9797958971132712*g1*
      g1p*Tr2U141*UNITMATRIX(3) + 2.065591117977289*g1*Tr31*UNITMATRIX(3) -
      3.794733192202055*g1p*Tr34*UNITMATRIX(3) + 53.333333333333336*Power(g3,4)
      *AbsSqr(MassG)*UNITMATRIX(3) + 0.5333333333333333*Tr2U111*Sqr(g1)*
      UNITMATRIX(3) + 1.8*Tr2U144*Sqr(g1p)*UNITMATRIX(3) + 2.8444444444444446*
      AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 1.4222222222222223*MassB*
      Conj(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 9.6*AbsSqr(MassG)*Sqr(g1p)*
      Sqr(g3)*UNITMATRIX(3) + 4.8*MassBp*Conj(MassG)*Sqr(g1p)*Sqr(g3)*
      UNITMATRIX(3) + 0.008888888888888889*Conj(MassB)*Sqr(g1)*(90*(2*MassB*(
      gD.transpose()*gD.conjugate()) - (TgD).transpose()*gD.conjugate()) + (
      1752*MassB*Sqr(g1) + 81*(2*MassB + MassBp)*Sqr(g1p) + 160*(2*MassB +
      MassG)*Sqr(g3))*UNITMATRIX(3)) + 0.01*Conj(MassBp)*Sqr(g1p)*(10*(-16*
      MassBp*(gD.transpose()*gD.conjugate()) + 2*MassBp*(-5 + Sqr(QS))*((Kappa)
      .transpose()*Kappa.conjugate()) + 8*((TgD).transpose()*gD.conjugate()) +
      5*((TKappa).transpose()*Kappa.conjugate()) - Sqr(QS)*((TKappa).transpose(
      )*Kappa.conjugate())) + 3*(24*(MassB + 2*MassBp)*Sqr(g1) + 160*(2*MassBp
      + MassG)*Sqr(g3) + 9*MassBp*Sqr(g1p)*(197 + Sqr(QS)))*UNITMATRIX(3))))
      .real();


   return beta_mDxbar2;
}

} // namespace flexiblesusy
