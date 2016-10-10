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

// File generated at Wed 3 Jun 2015 23:43:31

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of msbar2.
 *
 * @return one-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_msbar2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_msbar2;

   beta_msbar2 = Re(oneOver16PiSqr*(-0.31622776601683794*g1p*QS*Tr14 + 2*
      (mphi2 + ms2 + msbar2)*AbsSqr(Sigmax) + 2*AbsSqr(TSigmax) - 0.2*AbsSqr(
      MassBp)*Sqr(g1p)*Sqr(QS)));


   return beta_msbar2;
}

/**
 * Calculates the two-loop beta function of msbar2.
 *
 * @return two-loop beta function
 */
double CSE6SSM_soft_parameters::calc_beta_msbar2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
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
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_msbar2;

   beta_msbar2 = Re(twoLoop*(-1.2649110640673518*g1p*QS*Tr34 - 6*
      traceconjTKappaTpTKappa*AbsSqr(Sigmax) - 4*traceconjTLambda12TpTLambda12*
      AbsSqr(Sigmax) - 6*mphi2*traceKappaAdjKappa*AbsSqr(Sigmax) - 12*ms2*
      traceKappaAdjKappa*AbsSqr(Sigmax) - 6*msbar2*traceKappaAdjKappa*AbsSqr(
      Sigmax) - 6*traceKappaAdjKappaconjmDx2*AbsSqr(Sigmax) - 6*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Sigmax) - 4*mphi2*
      traceLambda12AdjLambda12*AbsSqr(Sigmax) - 8*ms2*traceLambda12AdjLambda12*
      AbsSqr(Sigmax) - 4*msbar2*traceLambda12AdjLambda12*AbsSqr(Sigmax) - 4*
      traceLambda12AdjLambda12conjmH2I2*AbsSqr(Sigmax) - 4*
      tracemH1I2AdjLambda12Lambda12*AbsSqr(Sigmax) - 4*mHd2*AbsSqr(Lambdax)*
      AbsSqr(Sigmax) - 4*mHu2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 4*mphi2*AbsSqr(
      Lambdax)*AbsSqr(Sigmax) - 8*ms2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 4*msbar2
      *AbsSqr(Lambdax)*AbsSqr(Sigmax) - 4*mHp2*AbsSqr(Sigmax)*AbsSqr(SigmaL) -
      4*mHpbar2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 8*mphi2*AbsSqr(Sigmax)*AbsSqr(
      SigmaL) - 4*ms2*AbsSqr(Sigmax)*AbsSqr(SigmaL) - 4*msbar2*AbsSqr(Sigmax)*
      AbsSqr(SigmaL) - 4*AbsSqr(Sigmax)*AbsSqr(TKappaPr) - 4*AbsSqr(Sigmax)*
      AbsSqr(TLambdax) - 6*traceKappaAdjKappa*AbsSqr(TSigmax) - 4*
      traceLambda12AdjLambda12*AbsSqr(TSigmax) - 4*AbsSqr(Lambdax)*AbsSqr(
      TSigmax) - 16*AbsSqr(Sigmax)*AbsSqr(TSigmax) - 4*AbsSqr(SigmaL)*AbsSqr(
      TSigmax) - 4*AbsSqr(Sigmax)*AbsSqr(TSigmaL) - 6*traceAdjKappaTKappa*Conj(
      TSigmax)*Sigmax - 4*traceAdjLambda12TLambda12*Conj(TSigmax)*Sigmax + 0.2*
      Tr2U144*Sqr(g1p)*Sqr(QS) + 0.06*Power(g1p,4)*AbsSqr(MassBp)*Sqr(QS)*(94 +
      Sqr(QS)) - 8*mphi2*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 8*ms2*Sqr(Conj(Sigmax
      ))*Sqr(Sigmax) - 8*msbar2*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 4*Conj(Lambdax)
      *Conj(TSigmax)*Sigmax*TLambdax - 6*traceconjTKappaTpKappa*Conj(Sigmax)*
      TSigmax - 4*traceconjTLambda12TpLambda12*Conj(Sigmax)*TSigmax - 4*Conj(
      Sigmax)*Conj(TKappaPr)*KappaPr*TSigmax - 4*Conj(Sigmax)*Conj(TLambdax)*
      Lambdax*TSigmax - 4*Conj(Sigmax)*Conj(TSigmaL)*SigmaL*TSigmax - 4*Conj(
      KappaPr)*((4*mphi2 + ms2 + msbar2)*AbsSqr(Sigmax)*KappaPr + Conj(TSigmax)
      *(Sigmax*TKappaPr + KappaPr*TSigmax)) - 4*Conj(SigmaL)*Conj(TSigmax)*
      Sigmax*TSigmaL));


   return beta_msbar2;
}

} // namespace flexiblesusy
