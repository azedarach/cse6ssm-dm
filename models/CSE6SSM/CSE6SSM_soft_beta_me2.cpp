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

// File generated at Wed 3 Jun 2015 23:43:28

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of me2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_me2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHp2*(hE*hE.adjoint()) + 4*mHd2*(Ye*
      Ye.adjoint()) + 4*(ThE*(ThE).adjoint()) + 4*(TYe*(TYe).adjoint()) + 4*(hE
      *mH1I2*hE.adjoint()) + 2*(hE*hE.adjoint()*me2) + 2*(me2*hE*hE.adjoint())
      + 2*(me2*Ye*Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*Ye.adjoint()*
      me2) + 1.5491933384829668*g1*Tr11*UNITMATRIX(3) + 0.31622776601683794*g1p
      *Tr14*UNITMATRIX(3) - 4.8*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 0.2*
      AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the two-loop beta function of me2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CSE6SSM_soft_parameters::calc_beta_me2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjfdTfd = TRACE_STRUCT.traceAdjfdTfd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTgDTpgD = TRACE_STRUCT.traceconjTgDTpgD;
   const double traceconjThETphE = TRACE_STRUCT.traceconjThETphE;
   const double traceconjTfdTpfd = TRACE_STRUCT.traceconjTfdTpfd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (twoLoop*(-12*traceconjTgDTpTgD*(hE*hE.adjoint()) - 4*
      traceconjThETpThE*(hE*hE.adjoint()) - 24*mHp2*tracegDAdjgD*(hE*hE.adjoint
      ()) - 12*tracegDAdjgDconjmq2*(hE*hE.adjoint()) - 12*
      tracegDconjmDxbar2AdjgD*(hE*hE.adjoint()) - 8*mHp2*tracehEAdjhE*(hE*
      hE.adjoint()) - 4*tracehEAdjhEme2*(hE*hE.adjoint()) - 4*tracehEmH1I2AdjhE
      *(hE*hE.adjoint()) - 8*mHp2*AbsSqr(SigmaL)*(hE*hE.adjoint()) - 4*mHpbar2*
      AbsSqr(SigmaL)*(hE*hE.adjoint()) - 4*mphi2*AbsSqr(SigmaL)*(hE*hE.adjoint(
      )) - 4*AbsSqr(TSigmaL)*(hE*hE.adjoint()) - 2.4*mHp2*Sqr(g1)*(hE*
      hE.adjoint()) + 2.4*mHp2*Sqr(g1p)*(hE*hE.adjoint()) + 12*mHp2*Sqr(g2)*(hE
      *hE.adjoint()) + 24*AbsSqr(MassWB)*Sqr(g2)*(hE*hE.adjoint()) - 12*
      traceAdjgDTgD*(hE*(ThE).adjoint()) - 4*traceAdjhEThE*(hE*(ThE).adjoint())
      + 2.4*MassB*Sqr(g1)*(hE*(ThE).adjoint()) - 2.4*MassBp*Sqr(g1p)*(hE*(ThE)
      .adjoint()) - 12*MassWB*Sqr(g2)*(hE*(ThE).adjoint()) - 4*Conj(SigmaL)*
      TSigmaL*(hE*(ThE).adjoint()) - 4*traceconjTfdTpTfd*(Ye*Ye.adjoint()) - 12
      *traceconjTYdTpTYd*(Ye*Ye.adjoint()) - 4*traceconjTYeTpTYe*(Ye*Ye.adjoint
      ()) - 8*mHd2*tracefdAdjfd*(Ye*Ye.adjoint()) - 4*tracefdAdjfdconjmSI2*(Ye*
      Ye.adjoint()) - 4*tracefdmH2I2Adjfd*(Ye*Ye.adjoint()) - 12*
      tracemd2YdAdjYd*(Ye*Ye.adjoint()) - 4*traceme2YeAdjYe*(Ye*Ye.adjoint()) -
      4*traceml2AdjYeYe*(Ye*Ye.adjoint()) - 12*tracemq2AdjYdYd*(Ye*Ye.adjoint(
      )) - 24*mHd2*traceYdAdjYd*(Ye*Ye.adjoint()) - 8*mHd2*traceYeAdjYe*(Ye*
      Ye.adjoint()) - 8*mHd2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()) - 4*mHu2*AbsSqr(
      Lambdax)*(Ye*Ye.adjoint()) - 4*ms2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()) - 4*
      AbsSqr(TLambdax)*(Ye*Ye.adjoint()) - 2.4*mHd2*Sqr(g1)*(Ye*Ye.adjoint()) +
      2.4*mHd2*Sqr(g1p)*(Ye*Ye.adjoint()) + 12*mHd2*Sqr(g2)*(Ye*Ye.adjoint())
      + 24*AbsSqr(MassWB)*Sqr(g2)*(Ye*Ye.adjoint()) - 4*traceAdjfdTfd*(Ye*(TYe)
      .adjoint()) - 12*traceAdjYdTYd*(Ye*(TYe).adjoint()) - 4*traceAdjYeTYe*(Ye
      *(TYe).adjoint()) + 2.4*MassB*Sqr(g1)*(Ye*(TYe).adjoint()) - 2.4*MassBp*
      Sqr(g1p)*(Ye*(TYe).adjoint()) - 12*MassWB*Sqr(g2)*(Ye*(TYe).adjoint()) -
      4*Conj(Lambdax)*TLambdax*(Ye*(TYe).adjoint()) - 12*traceconjTgDTpgD*(ThE*
      hE.adjoint()) - 4*traceconjThETphE*(ThE*hE.adjoint()) - 4*Conj(TSigmaL)*
      SigmaL*(ThE*hE.adjoint()) - 12*Conj(MassWB)*Sqr(g2)*(ThE*hE.adjoint()) -
      12*tracegDAdjgD*(ThE*(ThE).adjoint()) - 4*tracehEAdjhE*(ThE*(ThE).adjoint
      ()) - 4*AbsSqr(SigmaL)*(ThE*(ThE).adjoint()) - 2.4*Sqr(g1)*(ThE*(ThE)
      .adjoint()) + 2.4*Sqr(g1p)*(ThE*(ThE).adjoint()) + 12*Sqr(g2)*(ThE*(ThE)
      .adjoint()) - 4*traceconjTfdTpfd*(TYe*Ye.adjoint()) - 12*traceconjTYdTpYd
      *(TYe*Ye.adjoint()) - 4*traceconjTYeTpYe*(TYe*Ye.adjoint()) - 4*Conj(
      TLambdax)*Lambdax*(TYe*Ye.adjoint()) - 12*Conj(MassWB)*Sqr(g2)*(TYe*
      Ye.adjoint()) - 4*tracefdAdjfd*(TYe*(TYe).adjoint()) - 12*traceYdAdjYd*(
      TYe*(TYe).adjoint()) - 4*traceYeAdjYe*(TYe*(TYe).adjoint()) - 4*AbsSqr(
      Lambdax)*(TYe*(TYe).adjoint()) - 2.4*Sqr(g1)*(TYe*(TYe).adjoint()) + 2.4*
      Sqr(g1p)*(TYe*(TYe).adjoint()) + 12*Sqr(g2)*(TYe*(TYe).adjoint()) - 12*
      tracegDAdjgD*(hE*mH1I2*hE.adjoint()) - 4*tracehEAdjhE*(hE*mH1I2*
      hE.adjoint()) - 4*AbsSqr(SigmaL)*(hE*mH1I2*hE.adjoint()) - 2.4*Sqr(g1)*(
      hE*mH1I2*hE.adjoint()) + 2.4*Sqr(g1p)*(hE*mH1I2*hE.adjoint()) + 12*Sqr(g2
      )*(hE*mH1I2*hE.adjoint()) - 6*tracegDAdjgD*(hE*hE.adjoint()*me2) - 2*
      tracehEAdjhE*(hE*hE.adjoint()*me2) - 2*AbsSqr(SigmaL)*(hE*hE.adjoint()*
      me2) - 1.2*Sqr(g1)*(hE*hE.adjoint()*me2) + 1.2*Sqr(g1p)*(hE*hE.adjoint()*
      me2) + 6*Sqr(g2)*(hE*hE.adjoint()*me2) - 6*tracegDAdjgD*(me2*hE*
      hE.adjoint()) - 2*tracehEAdjhE*(me2*hE*hE.adjoint()) - 2*AbsSqr(SigmaL)*(
      me2*hE*hE.adjoint()) - 1.2*Sqr(g1)*(me2*hE*hE.adjoint()) + 1.2*Sqr(g1p)*(
      me2*hE*hE.adjoint()) + 6*Sqr(g2)*(me2*hE*hE.adjoint()) - 2*tracefdAdjfd*(
      me2*Ye*Ye.adjoint()) - 6*traceYdAdjYd*(me2*Ye*Ye.adjoint()) - 2*
      traceYeAdjYe*(me2*Ye*Ye.adjoint()) - 2*AbsSqr(Lambdax)*(me2*Ye*Ye.adjoint
      ()) - 1.2*Sqr(g1)*(me2*Ye*Ye.adjoint()) + 1.2*Sqr(g1p)*(me2*Ye*Ye.adjoint
      ()) + 6*Sqr(g2)*(me2*Ye*Ye.adjoint()) - 4*tracefdAdjfd*(Ye*ml2*Ye.adjoint
      ()) - 12*traceYdAdjYd*(Ye*ml2*Ye.adjoint()) - 4*traceYeAdjYe*(Ye*ml2*
      Ye.adjoint()) - 4*AbsSqr(Lambdax)*(Ye*ml2*Ye.adjoint()) - 2.4*Sqr(g1)*(Ye
      *ml2*Ye.adjoint()) + 2.4*Sqr(g1p)*(Ye*ml2*Ye.adjoint()) + 12*Sqr(g2)*(Ye*
      ml2*Ye.adjoint()) - 2*tracefdAdjfd*(Ye*Ye.adjoint()*me2) - 6*traceYdAdjYd
      *(Ye*Ye.adjoint()*me2) - 2*traceYeAdjYe*(Ye*Ye.adjoint()*me2) - 2*AbsSqr(
      Lambdax)*(Ye*Ye.adjoint()*me2) - 1.2*Sqr(g1)*(Ye*Ye.adjoint()*me2) + 1.2*
      Sqr(g1p)*(Ye*Ye.adjoint()*me2) + 6*Sqr(g2)*(Ye*Ye.adjoint()*me2) - 4*mHp2
      *(hE*fu.adjoint()*fu*hE.adjoint()) - 4*mHu2*(hE*fu.adjoint()*fu*
      hE.adjoint()) - 4*(hE*fu.adjoint()*Tfu*(ThE).adjoint()) - 8*mHp2*(hE*
      hE.adjoint()*hE*hE.adjoint()) - 4*(hE*hE.adjoint()*ThE*(ThE).adjoint()) -
      4*mHp2*(hE*(Lambda12).adjoint()*Lambda12*hE.adjoint()) - 4*ms2*(hE*(
      Lambda12).adjoint()*Lambda12*hE.adjoint()) - 4*(hE*(Lambda12).adjoint()*
      TLambda12*(ThE).adjoint()) - 4*(hE*(Tfu).adjoint()*Tfu*hE.adjoint()) - 4*
      (hE*(ThE).adjoint()*ThE*hE.adjoint()) - 4*(hE*(TLambda12).adjoint()*
      TLambda12*hE.adjoint()) - 8*mHd2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 4*(
      Ye*Ye.adjoint()*TYe*(TYe).adjoint()) - 4*(Ye*(TYe).adjoint()*TYe*
      Ye.adjoint()) - 4*(ThE*fu.adjoint()*fu*(ThE).adjoint()) - 4*(ThE*
      hE.adjoint()*hE*(ThE).adjoint()) - 4*(ThE*(Lambda12).adjoint()*Lambda12*(
      ThE).adjoint()) - 4*(ThE*(Tfu).adjoint()*fu*hE.adjoint()) - 4*(ThE*(ThE)
      .adjoint()*hE*hE.adjoint()) - 4*(ThE*(TLambda12).adjoint()*Lambda12*
      hE.adjoint()) - 4*(TYe*Ye.adjoint()*Ye*(TYe).adjoint()) - 4*(TYe*(TYe)
      .adjoint()*Ye*Ye.adjoint()) - 4*(hE*mH1I2*fu.adjoint()*fu*hE.adjoint()) -
      4*(hE*mH1I2*hE.adjoint()*hE*hE.adjoint()) - 4*(hE*mH1I2*(Lambda12)
      .adjoint()*Lambda12*hE.adjoint()) - 4*(hE*fu.adjoint()*fu*mH1I2*
      hE.adjoint()) - 2*(hE*fu.adjoint()*fu*hE.adjoint()*me2) - 4*(hE*
      fu.adjoint()*mSI2.conjugate()*fu*hE.adjoint()) - 4*(hE*hE.adjoint()*hE*
      mH1I2*hE.adjoint()) - 2*(hE*hE.adjoint()*hE*hE.adjoint()*me2) - 4*(hE*
      hE.adjoint()*me2*hE*hE.adjoint()) - 4*(hE*(Lambda12).adjoint()*
      mH2I2.conjugate()*Lambda12*hE.adjoint()) - 4*(hE*(Lambda12).adjoint()*
      Lambda12*mH1I2*hE.adjoint()) - 2*(hE*(Lambda12).adjoint()*Lambda12*
      hE.adjoint()*me2) - 2*(me2*hE*fu.adjoint()*fu*hE.adjoint()) - 2*(me2*hE*
      hE.adjoint()*hE*hE.adjoint()) - 2*(me2*hE*(Lambda12).adjoint()*Lambda12*
      hE.adjoint()) - 2*(me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 4*(Ye*ml2*
      Ye.adjoint()*Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()*me2*Ye*Ye.adjoint()) -
      4*(Ye*Ye.adjoint()*Ye*ml2*Ye.adjoint()) - 2*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*me2) + 0.9797958971132712*g1*g1p*Tr2U114*UNITMATRIX(3) +
      0.9797958971132712*g1*g1p*Tr2U141*UNITMATRIX(3) + 6.196773353931867*g1*
      Tr31*UNITMATRIX(3) + 1.2649110640673518*g1p*Tr34*UNITMATRIX(3) + 4.8*
      Tr2U111*Sqr(g1)*UNITMATRIX(3) + 0.2*Tr2U144*Sqr(g1p)*UNITMATRIX(3) + 0.24
      *Conj(MassB)*Sqr(g1)*(10*(-2*MassB*(hE*hE.adjoint()) - 2*MassB*(Ye*
      Ye.adjoint()) + ThE*hE.adjoint() + TYe*Ye.adjoint()) + (648*MassB*Sqr(g1)
      - (2*MassB + MassBp)*Sqr(g1p))*UNITMATRIX(3)) + 0.03*Conj(MassBp)*Sqr(
      g1p)*(80*(2*MassBp*(hE*hE.adjoint()) + 2*MassBp*(Ye*Ye.adjoint()) - ThE*
      hE.adjoint() - TYe*Ye.adjoint()) + (-8*(MassB + 2*MassBp)*Sqr(g1) +
      MassBp*Sqr(g1p)*(189 + Sqr(QS)))*UNITMATRIX(3)))).real();


   return beta_me2;
}

} // namespace flexiblesusy
