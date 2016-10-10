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

// File generated at Wed 3 Jun 2015 23:42:42

#include "CSE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME CSE6SSM_susy_parameters
#define TRACE_STRUCT susy_traces

CSE6SSM_susy_parameters::CSE6SSM_susy_parameters()
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), hE(Eigen::Matrix<double,3,2>::Zero(
   )), Ye(Eigen::Matrix<double,3,3>::Zero()), SigmaL(0), KappaPr(0), Sigmax(0),
   gD(Eigen::Matrix<double,3,3>::Zero()), Kappa(Eigen::Matrix<double,3,3>
   ::Zero()), Lambda12(Eigen::Matrix<double,2,2>::Zero()), Lambdax(0), fu(
   Eigen::Matrix<double,3,2>::Zero()), fd(Eigen::Matrix<double,3,2>::Zero()),
   Yu(Eigen::Matrix<double,3,3>::Zero()), MuPr(0), MuPhi(0), XiF(0), g1(0), g2(
   0), g3(0), g1p(0), vd(0), vu(0), vs(0), vsb(0), vphi(0), QS(0)
{
   set_number_of_parameters(numberOfParameters);
}

CSE6SSM_susy_parameters::CSE6SSM_susy_parameters(
   double scale_, double loops_, double thresholds_,
   const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,2>& hE_
   , const Eigen::Matrix<double,3,3>& Ye_, double SigmaL_, double KappaPr_,
   double Sigmax_, const Eigen::Matrix<double,3,3>& gD_, const Eigen::Matrix<
   double,3,3>& Kappa_, const Eigen::Matrix<double,2,2>& Lambda12_, double
   Lambdax_, const Eigen::Matrix<double,3,2>& fu_, const Eigen::Matrix<double,3
   ,2>& fd_, const Eigen::Matrix<double,3,3>& Yu_, double MuPr_, double MuPhi_,
   double XiF_, double g1_, double g2_, double g3_, double g1p_, double vd_,
   double vu_, double vs_, double vsb_, double vphi_, double QS_

)
   : Beta_function()
   , Yd(Yd_), hE(hE_), Ye(Ye_), SigmaL(SigmaL_), KappaPr(KappaPr_), Sigmax(
   Sigmax_), gD(gD_), Kappa(Kappa_), Lambda12(Lambda12_), Lambdax(Lambdax_), fu
   (fu_), fd(fd_), Yu(Yu_), MuPr(MuPr_), MuPhi(MuPhi_), XiF(XiF_), g1(g1_), g2(
   g2_), g3(g3_), g1p(g1p_), vd(vd_), vu(vu_), vs(vs_), vsb(vsb_), vphi(vphi_)
   , QS(QS_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd CSE6SSM_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

CSE6SSM_susy_parameters CSE6SSM_susy_parameters::calc_beta() const
{
   Susy_traces susy_traces;
   calc_susy_traces(susy_traces);

   Eigen::Matrix<double,3,3> beta_Yd(calc_beta_Yd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_hE(calc_beta_hE_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Ye(calc_beta_Ye_one_loop(TRACE_STRUCT));
   double beta_SigmaL(calc_beta_SigmaL_one_loop(TRACE_STRUCT));
   double beta_KappaPr(calc_beta_KappaPr_one_loop(TRACE_STRUCT));
   double beta_Sigmax(calc_beta_Sigmax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_gD(calc_beta_gD_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Kappa(calc_beta_Kappa_one_loop(TRACE_STRUCT))
      ;
   Eigen::Matrix<double,2,2> beta_Lambda12(calc_beta_Lambda12_one_loop(
      TRACE_STRUCT));
   double beta_Lambdax(calc_beta_Lambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_fu(calc_beta_fu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_fd(calc_beta_fd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yu(calc_beta_Yu_one_loop(TRACE_STRUCT));
   double beta_MuPr(calc_beta_MuPr_one_loop(TRACE_STRUCT));
   double beta_MuPhi(calc_beta_MuPhi_one_loop(TRACE_STRUCT));
   double beta_XiF(calc_beta_XiF_one_loop(TRACE_STRUCT));
   double beta_g1(calc_beta_g1_one_loop(TRACE_STRUCT));
   double beta_g2(calc_beta_g2_one_loop(TRACE_STRUCT));
   double beta_g3(calc_beta_g3_one_loop(TRACE_STRUCT));
   double beta_g1p(calc_beta_g1p_one_loop(TRACE_STRUCT));
   double beta_vd(calc_beta_vd_one_loop(TRACE_STRUCT));
   double beta_vu(calc_beta_vu_one_loop(TRACE_STRUCT));
   double beta_vs(calc_beta_vs_one_loop(TRACE_STRUCT));
   double beta_vsb(calc_beta_vsb_one_loop(TRACE_STRUCT));
   double beta_vphi(calc_beta_vphi_one_loop(TRACE_STRUCT));
   const double beta_QS = 0.;

   if (get_loops() > 1) {
      beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
      beta_hE += calc_beta_hE_two_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
      beta_SigmaL += calc_beta_SigmaL_two_loop(TRACE_STRUCT);
      beta_KappaPr += calc_beta_KappaPr_two_loop(TRACE_STRUCT);
      beta_Sigmax += calc_beta_Sigmax_two_loop(TRACE_STRUCT);
      beta_gD += calc_beta_gD_two_loop(TRACE_STRUCT);
      beta_Kappa += calc_beta_Kappa_two_loop(TRACE_STRUCT);
      beta_Lambda12 += calc_beta_Lambda12_two_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
      beta_fu += calc_beta_fu_two_loop(TRACE_STRUCT);
      beta_fd += calc_beta_fd_two_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
      beta_MuPr += calc_beta_MuPr_two_loop(TRACE_STRUCT);
      beta_MuPhi += calc_beta_MuPhi_two_loop(TRACE_STRUCT);
      beta_XiF += calc_beta_XiF_two_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
      beta_g1p += calc_beta_g1p_two_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_two_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_two_loop(TRACE_STRUCT);
      beta_vs += calc_beta_vs_two_loop(TRACE_STRUCT);
      beta_vsb += calc_beta_vsb_two_loop(TRACE_STRUCT);
      beta_vphi += calc_beta_vphi_two_loop(TRACE_STRUCT);

   }


   return CSE6SSM_susy_parameters(get_scale(), get_loops(), get_thresholds(),
                                  beta_Yd, beta_hE, beta_Ye, beta_SigmaL, beta_KappaPr, beta_Sigmax, beta_gD, beta_Kappa, beta_Lambda12, beta_Lambdax, beta_fu, beta_fd, beta_Yu, beta_MuPr, beta_MuPhi, beta_XiF, beta_g1, beta_g2, beta_g3, beta_g1p, beta_vd, beta_vu, beta_vs, beta_vsb, beta_vphi, beta_QS);
}

void CSE6SSM_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   hE = Eigen::Matrix<double,3,2>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   SigmaL = 0.;
   KappaPr = 0.;
   Sigmax = 0.;
   gD = Eigen::Matrix<double,3,3>::Zero();
   Kappa = Eigen::Matrix<double,3,3>::Zero();
   Lambda12 = Eigen::Matrix<double,2,2>::Zero();
   Lambdax = 0.;
   fu = Eigen::Matrix<double,3,2>::Zero();
   fd = Eigen::Matrix<double,3,2>::Zero();
   Yu = Eigen::Matrix<double,3,3>::Zero();
   MuPr = 0.;
   MuPhi = 0.;
   XiF = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   g1p = 0.;
   vd = 0.;
   vu = 0.;
   vs = 0.;
   vsb = 0.;
   vphi = 0.;
   QS = 0.;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu +
      gD.conjugate()*gD.transpose() - 0.016666666666666666*(2*Sqr(g1) + 3*Sqr(
      g1p) + 90*Sqr(g2) + 160*Sqr(g3))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(0.2*(-5*AbsSqr(SigmaL)*(gD.conjugate()*
         gD.transpose()) + 2*Sqr(g1)*(gD.conjugate()*gD.transpose()) + 3*Sqr(
         g1p)*(gD.conjugate()*gD.transpose()) - 10*(Yd.adjoint()*Yd*Yd.adjoint(
         )*Yd) - 10*(Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 10*(gD.conjugate()*
         gD.transpose()*gD.conjugate()*gD.transpose()) - 5*(gD.conjugate()*(
         Kappa).transpose()*Kappa.conjugate()*gD.transpose()) - 15*(
         gD.conjugate()*gD.transpose())*(gD*gD.adjoint()).trace() - 5*(
         gD.conjugate()*gD.transpose())*(hE*hE.adjoint()).trace() + Yd.adjoint(
         )*Yd*(-5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(g1p) - 5*(fd*fd.adjoint()
         ).trace() - 15*(Yd*Yd.adjoint()).trace() - 5*(Ye*Ye.adjoint()).trace()
         ) + Yu.adjoint()*Yu*(-5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(g1p) - 5*(fu
         *fu.adjoint()).trace() - 15*(Yu*Yu.adjoint()).trace())) +
         0.0002777777777777778*(1156*Power(g1,4) + 60*Sqr(g1p)*(9*Sqr(g2) + 16*
         Sqr(g3)) + Sqr(g1)*(-132*Sqr(g1p) + 360*Sqr(g2) + 640*Sqr(g3)) + 100*(
         297*Power(g2,4) + 256*Power(g3,4) + 288*Sqr(g2)*Sqr(g3)) + 9*Power(g1p
         ,4)*(189 + Sqr(QS)))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(Ye.adjoint()*Ye - 0.1*(3*Sqr(g1) + 2*Sqr(
      g1p) + 15*Sqr(g2))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (0.01*twoLoop*(10*(-20*(Ye.adjoint()*hE*hE.adjoint()*
         Ye + Ye.adjoint()*Ye*Ye.adjoint()*Ye) + Ye.adjoint()*Ye*(-10*AbsSqr(
         Lambdax) + 12*Sqr(g1) + 3*Sqr(g1p) - 10*(fd*fd.adjoint()).trace() - 30
         *(Yd*Yd.adjoint()).trace() - 10*(Ye*Ye.adjoint()).trace())) + (297*
         Power(g1,4) + 825*Power(g2,4) + 60*Sqr(g1p)*Sqr(g2) + 18*Sqr(g1)*(2*
         Sqr(g1p) + 5*Sqr(g2)) + Power(g1p,4)*(192 + Sqr(QS)))*UNITMATRIX(3)))
         .real();
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 0.45*Sqr(
      g1p) - 1.5*Sqr(g2) + (fd*fd.adjoint()).trace() + 3*(Yd*Yd.adjoint())
      .trace() + (Ye*Ye.adjoint()).trace()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.97*Power(g1,4) + 4.4325*Power(g1p,4) +
         8.25*Power(g2,4) - 0.09*Sqr(g1)*Sqr(g1p) + 0.9*Sqr(g1)*Sqr(g2) + 1.35*
         Sqr(g1p)*Sqr(g2) + 0.0225*Power(g1p,4)*Sqr(QS) - 3*Sqr(Conj(Lambdax))*
         Sqr(Lambdax) + Sqr(g1p)*(fd*fd.adjoint()).trace() - 0.4*Sqr(g1)*(Yd*
         Yd.adjoint()).trace() - 0.6*Sqr(g1p)*(Yd*Yd.adjoint()).trace() + 16*
         Sqr(g3)*(Yd*Yd.adjoint()).trace() + 1.2*Sqr(g1)*(Ye*Ye.adjoint())
         .trace() - 0.2*Sqr(g1p)*(Ye*Ye.adjoint()).trace() + 0.05*AbsSqr(
         Lambdax)*(-20*AbsSqr(Sigmax) - 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS) - 20*(fu*
         fu.adjoint()).trace() - 60*(Yu*Yu.adjoint()).trace() - 60*(Kappa*(
         Kappa).adjoint()).trace() - 40*(Lambda12*(Lambda12).adjoint()).trace()
         ) - 3*(fd*fd.adjoint()*fd*fd.adjoint()).trace() - 2*(fd*fd.adjoint()*
         fu*fu.adjoint()).trace() - 3*(gD*gD.adjoint()*Yd.transpose()*
         Yd.conjugate()).trace() - 2*(hE*hE.adjoint()*Ye*Ye.adjoint()).trace()
         - 9*(Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*
         Yd.adjoint()).trace() - 3*(Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace() -
         (Lambda12*(Lambda12).adjoint()*fd.transpose()*fd.conjugate()).trace())
         );
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 0.2*Sqr(
      g1p) - 1.5*Sqr(g2) + (fu*fu.adjoint()).trace() + 3*(Yu*Yu.adjoint())
      .trace()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.97*Power(g1,4) + 1.92*Power(g1p,4) +
         8.25*Power(g2,4) + 0.36*Sqr(g1)*Sqr(g1p) + 0.9*Sqr(g1)*Sqr(g2) + 0.6*
         Sqr(g1p)*Sqr(g2) + 0.01*Power(g1p,4)*Sqr(QS) - 3*Sqr(Conj(Lambdax))*
         Sqr(Lambdax) + 1.5*Sqr(g1p)*(fu*fu.adjoint()).trace() + 0.8*Sqr(g1)*(
         Yu*Yu.adjoint()).trace() - 0.3*Sqr(g1p)*(Yu*Yu.adjoint()).trace() + 16
         *Sqr(g3)*(Yu*Yu.adjoint()).trace() + 0.05*AbsSqr(Lambdax)*(-20*AbsSqr(
         Sigmax) + 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS) - 20*(fd*fd.adjoint()).trace()
         - 60*(Yd*Yd.adjoint()).trace() - 20*(Ye*Ye.adjoint()).trace() - 60*(
         Kappa*(Kappa).adjoint()).trace() - 40*(Lambda12*(Lambda12).adjoint())
         .trace()) - 2*(fd*fd.adjoint()*fu*fu.adjoint()).trace() - 3*(fu*
         fu.adjoint()*fu*fu.adjoint()).trace() - (fu*hE.adjoint()*hE*fu.adjoint
         ()).trace() - (fu*(Lambda12).adjoint()*Lambda12*fu.adjoint()).trace()
         - 3*(gD*gD.adjoint()*Yu.transpose()*Yu.conjugate()).trace() - 3*(Yd*
         Yu.adjoint()*Yu*Yd.adjoint()).trace() - 9*(Yu*Yu.adjoint()*Yu*
         Yu.adjoint()).trace()));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (0.06666666666666667*oneOver16PiSqr*(30*(Yd.conjugate()*
      Yd.transpose()) - (2*Sqr(g1) + 3*Sqr(g1p) + 40*Sqr(g3))*UNITMATRIX(3)))
      .real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yd.conjugate()*gD*gD.adjoint()*
         Yd.transpose() + Yd.conjugate()*Yd.transpose()*Yd.conjugate()*
         Yd.transpose() + Yd.conjugate()*Yu.transpose()*Yu.conjugate()*
         Yd.transpose()) + Yd.conjugate()*Yd.transpose()*(-2*AbsSqr(Lambdax) +
         0.4*Sqr(g1) + 0.6*Sqr(g1p) + 6*Sqr(g2) - 2*(fd*fd.adjoint()).trace() -
         6*(Yd*Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) +
         0.0011111111111111111*(1168*Power(g1,4) + 6400*Power(g3,4) + 960*Sqr(
         g1p)*Sqr(g3) + Sqr(g1)*(-96*Sqr(g1p) + 640*Sqr(g3)) + 9*Power(g1p,4)*(
         192 + Sqr(QS)))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.016666666666666666*(32*Sqr(g1) + 3*Sqr(g1p) + 160*Sqr(g3))*UNITMATRIX(3
      ))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-0.4*(5*(Yu.conjugate()*gD*gD.adjoint()*
         Yu.transpose() + Yu.conjugate()*Yd.transpose()*Yd.conjugate()*
         Yu.transpose() + Yu.conjugate()*Yu.transpose()*Yu.conjugate()*
         Yu.transpose()) + Yu.conjugate()*Yu.transpose()*(5*AbsSqr(Lambdax) +
         Sqr(g1) - Sqr(g1p) - 15*Sqr(g2) + 5*(fu*fu.adjoint()).trace() + 15*(Yu
         *Yu.adjoint()).trace())) + 0.0002777777777777778*(19456*Power(g1,4) +
         25600*Power(g3,4) + 960*Sqr(g1p)*Sqr(g3) + 256*Sqr(g1)*(3*Sqr(g1p) +
         40*Sqr(g3)) + 9*Power(g1p,4)*(189 + Sqr(QS)))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(hE.conjugate()*hE.transpose() +
      Ye.conjugate()*Ye.transpose()) - 0.05*(24*Sqr(g1) + Sqr(g1p))*UNITMATRIX(
      3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-0.4*(5*(hE.conjugate()*fu.transpose()*
         fu.conjugate()*hE.transpose() + hE.conjugate()*hE.transpose()*
         hE.conjugate()*hE.transpose() + hE.conjugate()*(Lambda12).transpose()*
         Lambda12.conjugate()*hE.transpose() + Ye.conjugate()*Ye.transpose()*
         Ye.conjugate()*Ye.transpose()) + hE.conjugate()*hE.transpose()*(5*
         AbsSqr(SigmaL) + 3*Sqr(g1) - 3*Sqr(g1p) - 15*Sqr(g2) + 15*(gD*
         gD.adjoint()).trace() + 5*(hE*hE.adjoint()).trace()) + Ye.conjugate()*
         Ye.transpose()*(5*AbsSqr(Lambdax) + 3*Sqr(g1) - 3*Sqr(g1p) - 15*Sqr(g2
         ) + 5*(fd*fd.adjoint()).trace() + 15*(Yd*Yd.adjoint()).trace() + 5*(Ye
         *Ye.adjoint()).trace())) + 0.0025*(5184*Power(g1,4) - 48*Sqr(g1)*Sqr(
         g1p) + Power(g1p,4)*(189 + Sqr(QS)))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(2*AbsSqr(Lambdax) + AbsSqr(Sigmax) - 0.05
      *Sqr(g1p)*Sqr(QS) + 3*(Kappa*(Kappa).adjoint()).trace() + 2*(Lambda12*(
      Lambda12).adjoint()).trace()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(0.005*Power(g1p,4)*Power(QS,4) - 2*AbsSqr
         (KappaPr)*AbsSqr(Sigmax) - 2*AbsSqr(Sigmax)*AbsSqr(SigmaL) + 0.47*
         Power(g1p,4)*Sqr(QS) - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 2*Sqr(Conj(
         Sigmax))*Sqr(Sigmax) + 0.1*AbsSqr(Lambdax)*(12*Sqr(g1) + 13*Sqr(g1p) +
         60*Sqr(g2) - Sqr(g1p)*Sqr(QS) - 20*(fd*fd.adjoint()).trace() - 20*(fu
         *fu.adjoint()).trace() - 60*(Yd*Yd.adjoint()).trace() - 20*(Ye*
         Ye.adjoint()).trace() - 60*(Yu*Yu.adjoint()).trace()) + 0.8*Sqr(g1)*(
         Kappa*(Kappa).adjoint()).trace() + 1.95*Sqr(g1p)*(Kappa*(Kappa)
         .adjoint()).trace() + 16*Sqr(g3)*(Kappa*(Kappa).adjoint()).trace() -
         0.15*Sqr(g1p)*Sqr(QS)*(Kappa*(Kappa).adjoint()).trace() + 1.2*Sqr(g1)*
         (Lambda12*(Lambda12).adjoint()).trace() + 1.3*Sqr(g1p)*(Lambda12*(
         Lambda12).adjoint()).trace() + 6*Sqr(g2)*(Lambda12*(Lambda12).adjoint(
         )).trace() - 0.1*Sqr(g1p)*Sqr(QS)*(Lambda12*(Lambda12).adjoint())
         .trace() - 2*(fu*(Lambda12).adjoint()*Lambda12*fu.adjoint()).trace() -
         6*(gD*(Kappa).adjoint()*Kappa*gD.adjoint()).trace() - 2*(hE*(Lambda12
         ).adjoint()*Lambda12*hE.adjoint()).trace() - 6*(Kappa*(Kappa).adjoint(
         )*Kappa*(Kappa).adjoint()).trace() - 4*(Lambda12*(Lambda12).adjoint()*
         Lambda12*(Lambda12).adjoint()).trace() - 2*(Lambda12*(Lambda12)
         .adjoint()*fd.transpose()*fd.conjugate()).trace()));
   }

   return anomDim;
}

double CLASSNAME::get_SsbarRSsbarR() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(Sigmax) - 0.05*Sqr(g1p)*Sqr(QS)));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(0.005*Power(g1p,4)*Power(QS,4) - 2*AbsSqr
         (KappaPr)*AbsSqr(Sigmax) - 2*AbsSqr(Lambdax)*AbsSqr(Sigmax) - 2*AbsSqr
         (Sigmax)*AbsSqr(SigmaL) + 0.47*Power(g1p,4)*Sqr(QS) - 2*Sqr(Conj(
         Sigmax))*Sqr(Sigmax) - 3*AbsSqr(Sigmax)*(Kappa*(Kappa).adjoint())
         .trace() - 2*AbsSqr(Sigmax)*(Lambda12*(Lambda12).adjoint()).trace()));
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH1ISH1I() const
{
   Eigen::Matrix<double,2,2> anomDim;

   anomDim = (oneOver16PiSqr*(fu.adjoint()*fu + hE.adjoint()*hE + (
      Lambda12).adjoint()*Lambda12 - 0.15*(2*Sqr(g1) + 3*Sqr(g1p) + 10*Sqr(g2))
      *UNITMATRIX(2))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-(AbsSqr(SigmaL)*(hE.adjoint()*hE)) + 1.2*
         Sqr(g1)*(hE.adjoint()*hE) - 0.2*Sqr(g1p)*(hE.adjoint()*hE) - 2*AbsSqr(
         Lambdax)*((Lambda12).adjoint()*Lambda12) - AbsSqr(Sigmax)*((Lambda12)
         .adjoint()*Lambda12) - 0.25*Sqr(g1p)*((Lambda12).adjoint()*Lambda12) +
         0.05*Sqr(g1p)*Sqr(QS)*((Lambda12).adjoint()*Lambda12) - 2*(fu.adjoint
         ()*fd*fd.adjoint()*fu) - 2*(fu.adjoint()*fu*fu.adjoint()*fu) - 2*(
         hE.adjoint()*hE*hE.adjoint()*hE) - 2*(hE.adjoint()*Ye*Ye.adjoint()*hE)
         - (Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12 - (
         Lambda12).adjoint()*fd.transpose()*fd.conjugate()*Lambda12 - 3*(
         hE.adjoint()*hE)*(gD*gD.adjoint()).trace() - hE.adjoint()*hE*(hE*
         hE.adjoint()).trace() + fu.adjoint()*fu*(-AbsSqr(Lambdax) + Sqr(g1p) -
         (fu*fu.adjoint()).trace() - 3*(Yu*Yu.adjoint()).trace()) - 3*((
         Lambda12).adjoint()*Lambda12)*(Kappa*(Kappa).adjoint()).trace() - 2*((
         Lambda12).adjoint()*Lambda12)*(Lambda12*(Lambda12).adjoint()).trace()
         + 0.0075*(396*Power(g1,4) + 1100*Power(g2,4) - 12*Sqr(g1)*(Sqr(g1p) -
         10*Sqr(g2)) + 180*Sqr(g1p)*Sqr(g2) + 3*Power(g1p,4)*(197 + Sqr(QS)))*
         UNITMATRIX(2))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH2ISH2I() const
{
   Eigen::Matrix<double,2,2> anomDim;

   anomDim = (oneOver16PiSqr*(fd.adjoint()*fd + Lambda12.conjugate()*(
      Lambda12).transpose() - 0.1*(3*Sqr(g1) + 2*Sqr(g1p) + 15*Sqr(g2))*
      UNITMATRIX(2))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*AbsSqr(Lambdax)*(Lambda12.conjugate()*(
         Lambda12).transpose()) - AbsSqr(Sigmax)*(Lambda12.conjugate()*(
         Lambda12).transpose()) + 0.25*Sqr(g1p)*(Lambda12.conjugate()*(Lambda12
         ).transpose()) + 0.05*Sqr(g1p)*Sqr(QS)*(Lambda12.conjugate()*(Lambda12
         ).transpose()) - 2*(fd.adjoint()*fd*fd.adjoint()*fd) - 2*(fd.adjoint()
         *fu*fu.adjoint()*fd) - Lambda12.conjugate()*fu.transpose()*
         fu.conjugate()*(Lambda12).transpose() - Lambda12.conjugate()*
         hE.transpose()*hE.conjugate()*(Lambda12).transpose() -
         Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(
         Lambda12).transpose() + fd.adjoint()*fd*(-AbsSqr(Lambdax) + 1.5*Sqr(
         g1p) - (fd*fd.adjoint()).trace() - 3*(Yd*Yd.adjoint()).trace() - (Ye*
         Ye.adjoint()).trace()) - 3*(Lambda12.conjugate()*(Lambda12).transpose(
         ))*(Kappa*(Kappa).adjoint()).trace() - 2*(Lambda12.conjugate()*(
         Lambda12).transpose())*(Lambda12*(Lambda12).adjoint()).trace() + 0.01*
         (297*Power(g1,4) + 825*Power(g2,4) + 60*Sqr(g1p)*Sqr(g2) + 18*Sqr(g1)*
         (2*Sqr(g1p) + 5*Sqr(g2)) + Power(g1p,4)*(192 + Sqr(QS)))*UNITMATRIX(2)
         )).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SSIRSSIR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(fd.conjugate()*fd.transpose() +
      fu.conjugate()*fu.transpose()) - 1.25*Sqr(g1p)*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(0.4*(-5*(fd.conjugate()*Lambda12*(Lambda12)
         .adjoint()*fd.transpose() + fd.conjugate()*fd.transpose()*fd.conjugate
         ()*fd.transpose() + fu.conjugate()*fu.transpose()*fu.conjugate()*
         fu.transpose() + fu.conjugate()*hE.transpose()*hE.conjugate()*
         fu.transpose() + fu.conjugate()*(Lambda12).transpose()*
         Lambda12.conjugate()*fu.transpose()) + fd.conjugate()*fd.transpose()*(
         -5*AbsSqr(Lambdax) + 3*Sqr(g1) - 3*Sqr(g1p) + 15*Sqr(g2) - 5*(fd*
         fd.adjoint()).trace() - 15*(Yd*Yd.adjoint()).trace() - 5*(Ye*
         Ye.adjoint()).trace()) + fu.conjugate()*fu.transpose()*(-5*AbsSqr(
         Lambdax) - 5*(fu*fu.adjoint()).trace() + 3*(Sqr(g1) - Sqr(g1p) + 5*Sqr
         (g2) - 5*(Yu*Yu.adjoint()).trace()))) + 0.0625*Power(g1p,4)*(213 + Sqr
         (QS))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxLSDxL() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(Kappa.conjugate()*(Kappa).transpose() -
      0.06666666666666667*(2*Sqr(g1) + 3*Sqr(g1p) + 40*Sqr(g3))*UNITMATRIX(3)))
      .real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Kappa.conjugate()*gD.transpose()*
         gD.conjugate()*(Kappa).transpose()) - Kappa.conjugate()*(Kappa)
         .transpose()*Kappa.conjugate()*(Kappa).transpose() + Kappa.conjugate()
         *(Kappa).transpose()*(-2*AbsSqr(Lambdax) - AbsSqr(Sigmax) + 0.25*Sqr(
         g1p) + 0.05*Sqr(g1p)*Sqr(QS) - 3*(Kappa*(Kappa).adjoint()).trace() - 2
         *(Lambda12*(Lambda12).adjoint()).trace()) + 0.0011111111111111111*(
         1168*Power(g1,4) + 6400*Power(g3,4) + 960*Sqr(g1p)*Sqr(g3) + Sqr(g1)*(
         -96*Sqr(g1p) + 640*Sqr(g3)) + 9*Power(g1p,4)*(192 + Sqr(QS)))*
         UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxbarRSDxbarR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(gD.adjoint()*gD) + (Kappa).adjoint()*
      Kappa - 0.016666666666666666*(8*Sqr(g1) + 27*Sqr(g1p) + 160*Sqr(g3))*
      UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*AbsSqr(Lambdax)*((Kappa).adjoint()*Kappa
         ) - AbsSqr(Sigmax)*((Kappa).adjoint()*Kappa) - 0.25*Sqr(g1p)*((Kappa)
         .adjoint()*Kappa) + 0.05*Sqr(g1p)*Sqr(QS)*((Kappa).adjoint()*Kappa) -
         2*(gD.adjoint()*gD*gD.adjoint()*gD) - 2*(gD.adjoint()*Yd.transpose()*
         Yd.conjugate()*gD) - 2*(gD.adjoint()*Yu.transpose()*Yu.conjugate()*gD)
         - (Kappa).adjoint()*Kappa*(Kappa).adjoint()*Kappa + 0.4*(gD.adjoint()
         *gD)*(-5*AbsSqr(SigmaL) + Sqr(g1) - Sqr(g1p) + 15*Sqr(g2) - 15*(gD*
         gD.adjoint()).trace() - 5*(hE*hE.adjoint()).trace()) - 3*((Kappa)
         .adjoint()*Kappa)*(Kappa*(Kappa).adjoint()).trace() - 2*((Kappa)
         .adjoint()*Kappa)*(Lambda12*(Lambda12).adjoint()).trace() +
         0.0002777777777777778*(4672*Power(g1,4) + 25600*Power(g3,4) + 8640*Sqr
         (g1p)*Sqr(g3) + 16*Sqr(g1)*(81*Sqr(g1p) + 160*Sqr(g3)) + 81*Power(g1p,
         4)*(197 + Sqr(QS)))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

double CLASSNAME::get_SHpSHp() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(SigmaL) - 0.3*Sqr(g1) - 0.2*Sqr(
      g1p) - 1.5*Sqr(g2) + 3*(gD*gD.adjoint()).trace() + (hE*hE.adjoint())
      .trace()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.97*Power(g1,4) + 1.92*Power(g1p,4) +
         8.25*Power(g2,4) - 2*AbsSqr(KappaPr)*AbsSqr(SigmaL) - AbsSqr(Sigmax)*
         AbsSqr(SigmaL) + 0.36*Sqr(g1)*Sqr(g1p) + 0.9*Sqr(g1)*Sqr(g2) + 0.6*Sqr
         (g1p)*Sqr(g2) + 0.01*Power(g1p,4)*Sqr(QS) - 3*Sqr(Conj(SigmaL))*Sqr(
         SigmaL) - 0.4*Sqr(g1)*(gD*gD.adjoint()).trace() + 0.9*Sqr(g1p)*(gD*
         gD.adjoint()).trace() + 16*Sqr(g3)*(gD*gD.adjoint()).trace() + 1.2*Sqr
         (g1)*(hE*hE.adjoint()).trace() + 0.3*Sqr(g1p)*(hE*hE.adjoint()).trace(
         ) - (fu*hE.adjoint()*hE*fu.adjoint()).trace() - 9*(gD*gD.adjoint()*gD*
         gD.adjoint()).trace() - 3*(gD*gD.adjoint()*Yd.transpose()*Yd.conjugate
         ()).trace() - 3*(gD*gD.adjoint()*Yu.transpose()*Yu.conjugate()).trace(
         ) - 3*(gD*(Kappa).adjoint()*Kappa*gD.adjoint()).trace() - 3*(hE*
         hE.adjoint()*hE*hE.adjoint()).trace() - 2*(hE*hE.adjoint()*Ye*
         Ye.adjoint()).trace() - (hE*(Lambda12).adjoint()*Lambda12*hE.adjoint()
         ).trace()));
   }

   return anomDim;
}

double CLASSNAME::get_SHpbarSHpbar() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(SigmaL) + 0.1*(-3*Sqr(g1) - 2*Sqr(
      g1p) - 15*Sqr(g2))));

   if (get_loops() > 1) {
      anomDim += Re(0.01*twoLoop*(297*Power(g1,4) + 192*Power(g1p,4) +
         825*Power(g2,4) - 200*AbsSqr(KappaPr)*AbsSqr(SigmaL) - 100*AbsSqr(
         Sigmax)*AbsSqr(SigmaL) + 36*Sqr(g1)*Sqr(g1p) + 90*Sqr(g1)*Sqr(g2) + 60
         *Sqr(g1p)*Sqr(g2) + Power(g1p,4)*Sqr(QS) - 300*Sqr(Conj(SigmaL))*Sqr(
         SigmaL) - 300*AbsSqr(SigmaL)*(gD*gD.adjoint()).trace() - 100*AbsSqr(
         SigmaL)*(hE*hE.adjoint()).trace()));
   }

   return anomDim;
}

double CLASSNAME::get_SphiRSphiR() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax) + 2*
      AbsSqr(SigmaL)));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(-4*AbsSqr(KappaPr)*(AbsSqr(Sigmax) + 2*
         AbsSqr(SigmaL)) - 8*Sqr(Conj(KappaPr))*Sqr(KappaPr) - 2*Sqr(Conj(
         Sigmax))*Sqr(Sigmax) + 0.4*AbsSqr(SigmaL)*(-10*AbsSqr(SigmaL) + 3*Sqr(
         g1) + 2*Sqr(g1p) + 15*Sqr(g2) - 15*(gD*gD.adjoint()).trace() - 5*(hE*
         hE.adjoint()).trace()) + Conj(Sigmax)*(-2*AbsSqr(Lambdax)*Sigmax + 0.1
         *Sigmax*Sqr(g1p)*Sqr(QS) - 3*Sigmax*(Kappa*(Kappa).adjoint()).trace()
         - 2*Sigmax*(Lambda12*(Lambda12).adjoint()).trace())));
   }

   return anomDim;
}


const Eigen::ArrayXd CSE6SSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = hE(0,0);
   pars(10) = hE(0,1);
   pars(11) = hE(1,0);
   pars(12) = hE(1,1);
   pars(13) = hE(2,0);
   pars(14) = hE(2,1);
   pars(15) = Ye(0,0);
   pars(16) = Ye(0,1);
   pars(17) = Ye(0,2);
   pars(18) = Ye(1,0);
   pars(19) = Ye(1,1);
   pars(20) = Ye(1,2);
   pars(21) = Ye(2,0);
   pars(22) = Ye(2,1);
   pars(23) = Ye(2,2);
   pars(24) = SigmaL;
   pars(25) = KappaPr;
   pars(26) = Sigmax;
   pars(27) = gD(0,0);
   pars(28) = gD(0,1);
   pars(29) = gD(0,2);
   pars(30) = gD(1,0);
   pars(31) = gD(1,1);
   pars(32) = gD(1,2);
   pars(33) = gD(2,0);
   pars(34) = gD(2,1);
   pars(35) = gD(2,2);
   pars(36) = Kappa(0,0);
   pars(37) = Kappa(0,1);
   pars(38) = Kappa(0,2);
   pars(39) = Kappa(1,0);
   pars(40) = Kappa(1,1);
   pars(41) = Kappa(1,2);
   pars(42) = Kappa(2,0);
   pars(43) = Kappa(2,1);
   pars(44) = Kappa(2,2);
   pars(45) = Lambda12(0,0);
   pars(46) = Lambda12(0,1);
   pars(47) = Lambda12(1,0);
   pars(48) = Lambda12(1,1);
   pars(49) = Lambdax;
   pars(50) = fu(0,0);
   pars(51) = fu(0,1);
   pars(52) = fu(1,0);
   pars(53) = fu(1,1);
   pars(54) = fu(2,0);
   pars(55) = fu(2,1);
   pars(56) = fd(0,0);
   pars(57) = fd(0,1);
   pars(58) = fd(1,0);
   pars(59) = fd(1,1);
   pars(60) = fd(2,0);
   pars(61) = fd(2,1);
   pars(62) = Yu(0,0);
   pars(63) = Yu(0,1);
   pars(64) = Yu(0,2);
   pars(65) = Yu(1,0);
   pars(66) = Yu(1,1);
   pars(67) = Yu(1,2);
   pars(68) = Yu(2,0);
   pars(69) = Yu(2,1);
   pars(70) = Yu(2,2);
   pars(71) = MuPr;
   pars(72) = MuPhi;
   pars(73) = XiF;
   pars(74) = g1;
   pars(75) = g2;
   pars(76) = g3;
   pars(77) = g1p;
   pars(78) = vd;
   pars(79) = vu;
   pars(80) = vs;
   pars(81) = vsb;
   pars(82) = vphi;
   pars(83) = QS;


   return pars;
}

void CSE6SSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "susy parameters:\n"
           "----------------------------------------\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "hE = " << hE << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "SigmaL = " << SigmaL << '\n';
   ostr << "KappaPr = " << KappaPr << '\n';
   ostr << "Sigmax = " << Sigmax << '\n';
   ostr << "gD = " << gD << '\n';
   ostr << "Kappa = " << Kappa << '\n';
   ostr << "Lambda12 = " << Lambda12 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "fu = " << fu << '\n';
   ostr << "fd = " << fd << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "MuPr = " << MuPr << '\n';
   ostr << "MuPhi = " << MuPhi << '\n';
   ostr << "XiF = " << XiF << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "g1p = " << g1p << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vs = " << vs << '\n';
   ostr << "vsb = " << vsb << '\n';
   ostr << "vphi = " << vphi << '\n';
   ostr << "QS = " << QS << '\n';

}

void CSE6SSM_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   Yd(0,0) = pars(0);
   Yd(0,1) = pars(1);
   Yd(0,2) = pars(2);
   Yd(1,0) = pars(3);
   Yd(1,1) = pars(4);
   Yd(1,2) = pars(5);
   Yd(2,0) = pars(6);
   Yd(2,1) = pars(7);
   Yd(2,2) = pars(8);
   hE(0,0) = pars(9);
   hE(0,1) = pars(10);
   hE(1,0) = pars(11);
   hE(1,1) = pars(12);
   hE(2,0) = pars(13);
   hE(2,1) = pars(14);
   Ye(0,0) = pars(15);
   Ye(0,1) = pars(16);
   Ye(0,2) = pars(17);
   Ye(1,0) = pars(18);
   Ye(1,1) = pars(19);
   Ye(1,2) = pars(20);
   Ye(2,0) = pars(21);
   Ye(2,1) = pars(22);
   Ye(2,2) = pars(23);
   SigmaL = pars(24);
   KappaPr = pars(25);
   Sigmax = pars(26);
   gD(0,0) = pars(27);
   gD(0,1) = pars(28);
   gD(0,2) = pars(29);
   gD(1,0) = pars(30);
   gD(1,1) = pars(31);
   gD(1,2) = pars(32);
   gD(2,0) = pars(33);
   gD(2,1) = pars(34);
   gD(2,2) = pars(35);
   Kappa(0,0) = pars(36);
   Kappa(0,1) = pars(37);
   Kappa(0,2) = pars(38);
   Kappa(1,0) = pars(39);
   Kappa(1,1) = pars(40);
   Kappa(1,2) = pars(41);
   Kappa(2,0) = pars(42);
   Kappa(2,1) = pars(43);
   Kappa(2,2) = pars(44);
   Lambda12(0,0) = pars(45);
   Lambda12(0,1) = pars(46);
   Lambda12(1,0) = pars(47);
   Lambda12(1,1) = pars(48);
   Lambdax = pars(49);
   fu(0,0) = pars(50);
   fu(0,1) = pars(51);
   fu(1,0) = pars(52);
   fu(1,1) = pars(53);
   fu(2,0) = pars(54);
   fu(2,1) = pars(55);
   fd(0,0) = pars(56);
   fd(0,1) = pars(57);
   fd(1,0) = pars(58);
   fd(1,1) = pars(59);
   fd(2,0) = pars(60);
   fd(2,1) = pars(61);
   Yu(0,0) = pars(62);
   Yu(0,1) = pars(63);
   Yu(0,2) = pars(64);
   Yu(1,0) = pars(65);
   Yu(1,1) = pars(66);
   Yu(1,2) = pars(67);
   Yu(2,0) = pars(68);
   Yu(2,1) = pars(69);
   Yu(2,2) = pars(70);
   MuPr = pars(71);
   MuPhi = pars(72);
   XiF = pars(73);
   g1 = pars(74);
   g2 = pars(75);
   g3 = pars(76);
   g1p = pars(77);
   vd = pars(78);
   vu = pars(79);
   vs = pars(80);
   vsb = pars(81);
   vphi = pars(82);
   QS = pars(83);

}

void CSE6SSM_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   TRACE_STRUCT.tracefdAdjfd = Re((fd*fd.adjoint()).trace());
   TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
   TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfu = Re((fu*fu.adjoint()).trace());
   TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappa = Re((Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12 = Re((Lambda12*(Lambda12).adjoint())
      .trace());
   TRACE_STRUCT.tracefdAdjfdfdAdjfd = Re((fd*fd.adjoint()*fd*fd.adjoint())
      .trace());
   TRACE_STRUCT.tracefdAdjfdfuAdjfu = Re((fd*fd.adjoint()*fu*fu.adjoint())
      .trace());
   TRACE_STRUCT.tracegDAdjgDTpYdconjYd = Re((gD*gD.adjoint()*Yd.transpose()*
      Yd.conjugate()).trace());
   TRACE_STRUCT.tracehEAdjhEYeAdjYe = Re((hE*hE.adjoint()*Ye*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceLambda12AdjLambda12Tpfdconjfd = Re((Lambda12*(Lambda12)
      .adjoint()*fd.transpose()*fd.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjgD = Re((gD*gD.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhE = Re((hE*hE.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjhEhEAdjfu = Re((fu*hE.adjoint()*hE*fu.adjoint())
      .trace());
   TRACE_STRUCT.tracegDAdjgDgDAdjgD = Re((gD*gD.adjoint()*gD*gD.adjoint())
      .trace());
   TRACE_STRUCT.tracegDAdjgDTpYuconjYu = Re((gD*gD.adjoint()*Yu.transpose()*
      Yu.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjKappaKappaAdjgD = Re((gD*(Kappa).adjoint()*Kappa*
      gD.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhEhEAdjhE = Re((hE*hE.adjoint()*hE*hE.adjoint())
      .trace());
   TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhE = Re((hE*(Lambda12).adjoint()*
      Lambda12*hE.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjLambda12Lambda12Adjfu = Re((fu*(Lambda12).adjoint()*
      Lambda12*fu.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = Re((Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12 = Re((Lambda12*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfufuAdjfu = Re((fu*fu.adjoint()*fu*fu.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
      .trace());

}

std::ostream& operator<<(std::ostream& ostr, const CSE6SSM_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
