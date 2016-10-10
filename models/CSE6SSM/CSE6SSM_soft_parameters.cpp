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

// File generated at Wed 3 Jun 2015 23:42:56

#include "CSE6SSM_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define TRACE_STRUCT soft_traces

CSE6SSM_soft_parameters::CSE6SSM_soft_parameters()
   : CSE6SSM_susy_parameters()
   , TYd(Eigen::Matrix<double,3,3>::Zero()), ThE(Eigen::Matrix<double,3,2>
   ::Zero()), TYe(Eigen::Matrix<double,3,3>::Zero()), TSigmaL(0), TKappaPr(0),
   TSigmax(0), TgD(Eigen::Matrix<double,3,3>::Zero()), TKappa(Eigen::Matrix<
   double,3,3>::Zero()), TLambda12(Eigen::Matrix<double,2,2>::Zero()), TLambdax
   (0), Tfu(Eigen::Matrix<double,3,2>::Zero()), Tfd(Eigen::Matrix<double,3,2>
   ::Zero()), TYu(Eigen::Matrix<double,3,3>::Zero()), BMuPr(0), BMuPhi(0), LXiF
   (0), mq2(Eigen::Matrix<double,3,3>::Zero()), ml2(Eigen::Matrix<double,3,3>
   ::Zero()), mHd2(0), mHu2(0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(
   Eigen::Matrix<double,3,3>::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()),
   ms2(0), msbar2(0), mH1I2(Eigen::Matrix<double,2,2>::Zero()), mH2I2(
   Eigen::Matrix<double,2,2>::Zero()), mSI2(Eigen::Matrix<double,3,3>::Zero()),
   mDx2(Eigen::Matrix<double,3,3>::Zero()), mDxbar2(Eigen::Matrix<double,3,3>
   ::Zero()), mHp2(0), mHpbar2(0), mphi2(0), MassB(0), MassWB(0), MassG(0),
   MassBp(0)

{
   set_number_of_parameters(numberOfParameters);
}

CSE6SSM_soft_parameters::CSE6SSM_soft_parameters(
   const CSE6SSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,2>&
   ThE_, const Eigen::Matrix<double,3,3>& TYe_, double TSigmaL_, double
   TKappaPr_, double TSigmax_, const Eigen::Matrix<double,3,3>& TgD_, const
   Eigen::Matrix<double,3,3>& TKappa_, const Eigen::Matrix<double,2,2>&
   TLambda12_, double TLambdax_, const Eigen::Matrix<double,3,2>& Tfu_, const
   Eigen::Matrix<double,3,2>& Tfd_, const Eigen::Matrix<double,3,3>& TYu_,
   double BMuPr_, double BMuPhi_, double LXiF_, const Eigen::Matrix<double,3,3>
   & mq2_, const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_,
   const Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_
   , const Eigen::Matrix<double,3,3>& me2_, double ms2_, double msbar2_, const
   Eigen::Matrix<double,2,2>& mH1I2_, const Eigen::Matrix<double,2,2>& mH2I2_,
   const Eigen::Matrix<double,3,3>& mSI2_, const Eigen::Matrix<double,3,3>&
   mDx2_, const Eigen::Matrix<double,3,3>& mDxbar2_, double mHp2_, double
   mHpbar2_, double mphi2_, double MassB_, double MassWB_, double MassG_,
   double MassBp_

)
   : CSE6SSM_susy_parameters(susy_model)
   , TYd(TYd_), ThE(ThE_), TYe(TYe_), TSigmaL(TSigmaL_), TKappaPr(TKappaPr_),
   TSigmax(TSigmax_), TgD(TgD_), TKappa(TKappa_), TLambda12(TLambda12_),
   TLambdax(TLambdax_), Tfu(Tfu_), Tfd(Tfd_), TYu(TYu_), BMuPr(BMuPr_), BMuPhi(
   BMuPhi_), LXiF(LXiF_), mq2(mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(
   md2_), mu2(mu2_), me2(me2_), ms2(ms2_), msbar2(msbar2_), mH1I2(mH1I2_),
   mH2I2(mH2I2_), mSI2(mSI2_), mDx2(mDx2_), mDxbar2(mDxbar2_), mHp2(mHp2_),
   mHpbar2(mHpbar2_), mphi2(mphi2_), MassB(MassB_), MassWB(MassWB_), MassG(
   MassG_), MassBp(MassBp_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd CSE6SSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

CSE6SSM_soft_parameters CSE6SSM_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   Eigen::Matrix<double,3,3> beta_TYd(calc_beta_TYd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_ThE(calc_beta_ThE_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYe(calc_beta_TYe_one_loop(TRACE_STRUCT));
   double beta_TSigmaL(calc_beta_TSigmaL_one_loop(TRACE_STRUCT));
   double beta_TKappaPr(calc_beta_TKappaPr_one_loop(TRACE_STRUCT));
   double beta_TSigmax(calc_beta_TSigmax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TgD(calc_beta_TgD_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TKappa(calc_beta_TKappa_one_loop(TRACE_STRUCT
      ));
   Eigen::Matrix<double,2,2> beta_TLambda12(calc_beta_TLambda12_one_loop(
      TRACE_STRUCT));
   double beta_TLambdax(calc_beta_TLambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_Tfu(calc_beta_Tfu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,2> beta_Tfd(calc_beta_Tfd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYu(calc_beta_TYu_one_loop(TRACE_STRUCT));
   double beta_BMuPr(calc_beta_BMuPr_one_loop(TRACE_STRUCT));
   double beta_BMuPhi(calc_beta_BMuPhi_one_loop(TRACE_STRUCT));
   double beta_LXiF(calc_beta_LXiF_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mq2(calc_beta_mq2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_ml2(calc_beta_ml2_one_loop(TRACE_STRUCT));
   double beta_mHd2(calc_beta_mHd2_one_loop(TRACE_STRUCT));
   double beta_mHu2(calc_beta_mHu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_md2(calc_beta_md2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_me2(calc_beta_me2_one_loop(TRACE_STRUCT));
   double beta_ms2(calc_beta_ms2_one_loop(TRACE_STRUCT));
   double beta_msbar2(calc_beta_msbar2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,2,2> beta_mH1I2(calc_beta_mH1I2_one_loop(TRACE_STRUCT))
      ;
   Eigen::Matrix<double,2,2> beta_mH2I2(calc_beta_mH2I2_one_loop(TRACE_STRUCT))
      ;
   Eigen::Matrix<double,3,3> beta_mSI2(calc_beta_mSI2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mDx2(calc_beta_mDx2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mDxbar2(calc_beta_mDxbar2_one_loop(
      TRACE_STRUCT));
   double beta_mHp2(calc_beta_mHp2_one_loop(TRACE_STRUCT));
   double beta_mHpbar2(calc_beta_mHpbar2_one_loop(TRACE_STRUCT));
   double beta_mphi2(calc_beta_mphi2_one_loop(TRACE_STRUCT));
   double beta_MassB(calc_beta_MassB_one_loop(TRACE_STRUCT));
   double beta_MassWB(calc_beta_MassWB_one_loop(TRACE_STRUCT));
   double beta_MassG(calc_beta_MassG_one_loop(TRACE_STRUCT));
   double beta_MassBp(calc_beta_MassBp_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_TYd += calc_beta_TYd_two_loop(TRACE_STRUCT);
      beta_ThE += calc_beta_ThE_two_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_two_loop(TRACE_STRUCT);
      beta_TSigmaL += calc_beta_TSigmaL_two_loop(TRACE_STRUCT);
      beta_TKappaPr += calc_beta_TKappaPr_two_loop(TRACE_STRUCT);
      beta_TSigmax += calc_beta_TSigmax_two_loop(TRACE_STRUCT);
      beta_TgD += calc_beta_TgD_two_loop(TRACE_STRUCT);
      beta_TKappa += calc_beta_TKappa_two_loop(TRACE_STRUCT);
      beta_TLambda12 += calc_beta_TLambda12_two_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_two_loop(TRACE_STRUCT);
      beta_Tfu += calc_beta_Tfu_two_loop(TRACE_STRUCT);
      beta_Tfd += calc_beta_Tfd_two_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_two_loop(TRACE_STRUCT);
      beta_BMuPr += calc_beta_BMuPr_two_loop(TRACE_STRUCT);
      beta_BMuPhi += calc_beta_BMuPhi_two_loop(TRACE_STRUCT);
      beta_LXiF += calc_beta_LXiF_two_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_two_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_two_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_two_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_two_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_two_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_two_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_two_loop(TRACE_STRUCT);
      beta_msbar2 += calc_beta_msbar2_two_loop(TRACE_STRUCT);
      beta_mH1I2 += calc_beta_mH1I2_two_loop(TRACE_STRUCT);
      beta_mH2I2 += calc_beta_mH2I2_two_loop(TRACE_STRUCT);
      beta_mSI2 += calc_beta_mSI2_two_loop(TRACE_STRUCT);
      beta_mDx2 += calc_beta_mDx2_two_loop(TRACE_STRUCT);
      beta_mDxbar2 += calc_beta_mDxbar2_two_loop(TRACE_STRUCT);
      beta_mHp2 += calc_beta_mHp2_two_loop(TRACE_STRUCT);
      beta_mHpbar2 += calc_beta_mHpbar2_two_loop(TRACE_STRUCT);
      beta_mphi2 += calc_beta_mphi2_two_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);
      beta_MassBp += calc_beta_MassBp_two_loop(TRACE_STRUCT);

   }


   const CSE6SSM_susy_parameters susy_betas(CSE6SSM_susy_parameters::calc_beta());

   return CSE6SSM_soft_parameters(susy_betas, beta_TYd, beta_ThE, beta_TYe, beta_TSigmaL, beta_TKappaPr, beta_TSigmax, beta_TgD, beta_TKappa, beta_TLambda12, beta_TLambdax, beta_Tfu, beta_Tfd, beta_TYu, beta_BMuPr, beta_BMuPhi, beta_LXiF, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_msbar2, beta_mH1I2, beta_mH2I2, beta_mSI2, beta_mDx2, beta_mDxbar2, beta_mHp2, beta_mHpbar2, beta_mphi2, beta_MassB, beta_MassWB, beta_MassG, beta_MassBp);
}

void CSE6SSM_soft_parameters::clear()
{
   CSE6SSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   ThE = Eigen::Matrix<double,3,2>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TSigmaL = 0.;
   TKappaPr = 0.;
   TSigmax = 0.;
   TgD = Eigen::Matrix<double,3,3>::Zero();
   TKappa = Eigen::Matrix<double,3,3>::Zero();
   TLambda12 = Eigen::Matrix<double,2,2>::Zero();
   TLambdax = 0.;
   Tfu = Eigen::Matrix<double,3,2>::Zero();
   Tfd = Eigen::Matrix<double,3,2>::Zero();
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMuPr = 0.;
   BMuPhi = 0.;
   LXiF = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   msbar2 = 0.;
   mH1I2 = Eigen::Matrix<double,2,2>::Zero();
   mH2I2 = Eigen::Matrix<double,2,2>::Zero();
   mSI2 = Eigen::Matrix<double,3,3>::Zero();
   mDx2 = Eigen::Matrix<double,3,3>::Zero();
   mDxbar2 = Eigen::Matrix<double,3,3>::Zero();
   mHp2 = 0.;
   mHpbar2 = 0.;
   mphi2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
   MassBp = 0.;

}

const Eigen::ArrayXd CSE6SSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(CSE6SSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(84) = TYd(0,0);
   pars(85) = TYd(0,1);
   pars(86) = TYd(0,2);
   pars(87) = TYd(1,0);
   pars(88) = TYd(1,1);
   pars(89) = TYd(1,2);
   pars(90) = TYd(2,0);
   pars(91) = TYd(2,1);
   pars(92) = TYd(2,2);
   pars(93) = ThE(0,0);
   pars(94) = ThE(0,1);
   pars(95) = ThE(1,0);
   pars(96) = ThE(1,1);
   pars(97) = ThE(2,0);
   pars(98) = ThE(2,1);
   pars(99) = TYe(0,0);
   pars(100) = TYe(0,1);
   pars(101) = TYe(0,2);
   pars(102) = TYe(1,0);
   pars(103) = TYe(1,1);
   pars(104) = TYe(1,2);
   pars(105) = TYe(2,0);
   pars(106) = TYe(2,1);
   pars(107) = TYe(2,2);
   pars(108) = TSigmaL;
   pars(109) = TKappaPr;
   pars(110) = TSigmax;
   pars(111) = TgD(0,0);
   pars(112) = TgD(0,1);
   pars(113) = TgD(0,2);
   pars(114) = TgD(1,0);
   pars(115) = TgD(1,1);
   pars(116) = TgD(1,2);
   pars(117) = TgD(2,0);
   pars(118) = TgD(2,1);
   pars(119) = TgD(2,2);
   pars(120) = TKappa(0,0);
   pars(121) = TKappa(0,1);
   pars(122) = TKappa(0,2);
   pars(123) = TKappa(1,0);
   pars(124) = TKappa(1,1);
   pars(125) = TKappa(1,2);
   pars(126) = TKappa(2,0);
   pars(127) = TKappa(2,1);
   pars(128) = TKappa(2,2);
   pars(129) = TLambda12(0,0);
   pars(130) = TLambda12(0,1);
   pars(131) = TLambda12(1,0);
   pars(132) = TLambda12(1,1);
   pars(133) = TLambdax;
   pars(134) = Tfu(0,0);
   pars(135) = Tfu(0,1);
   pars(136) = Tfu(1,0);
   pars(137) = Tfu(1,1);
   pars(138) = Tfu(2,0);
   pars(139) = Tfu(2,1);
   pars(140) = Tfd(0,0);
   pars(141) = Tfd(0,1);
   pars(142) = Tfd(1,0);
   pars(143) = Tfd(1,1);
   pars(144) = Tfd(2,0);
   pars(145) = Tfd(2,1);
   pars(146) = TYu(0,0);
   pars(147) = TYu(0,1);
   pars(148) = TYu(0,2);
   pars(149) = TYu(1,0);
   pars(150) = TYu(1,1);
   pars(151) = TYu(1,2);
   pars(152) = TYu(2,0);
   pars(153) = TYu(2,1);
   pars(154) = TYu(2,2);
   pars(155) = BMuPr;
   pars(156) = BMuPhi;
   pars(157) = LXiF;
   pars(158) = mq2(0,0);
   pars(159) = mq2(0,1);
   pars(160) = mq2(0,2);
   pars(161) = mq2(1,0);
   pars(162) = mq2(1,1);
   pars(163) = mq2(1,2);
   pars(164) = mq2(2,0);
   pars(165) = mq2(2,1);
   pars(166) = mq2(2,2);
   pars(167) = ml2(0,0);
   pars(168) = ml2(0,1);
   pars(169) = ml2(0,2);
   pars(170) = ml2(1,0);
   pars(171) = ml2(1,1);
   pars(172) = ml2(1,2);
   pars(173) = ml2(2,0);
   pars(174) = ml2(2,1);
   pars(175) = ml2(2,2);
   pars(176) = mHd2;
   pars(177) = mHu2;
   pars(178) = md2(0,0);
   pars(179) = md2(0,1);
   pars(180) = md2(0,2);
   pars(181) = md2(1,0);
   pars(182) = md2(1,1);
   pars(183) = md2(1,2);
   pars(184) = md2(2,0);
   pars(185) = md2(2,1);
   pars(186) = md2(2,2);
   pars(187) = mu2(0,0);
   pars(188) = mu2(0,1);
   pars(189) = mu2(0,2);
   pars(190) = mu2(1,0);
   pars(191) = mu2(1,1);
   pars(192) = mu2(1,2);
   pars(193) = mu2(2,0);
   pars(194) = mu2(2,1);
   pars(195) = mu2(2,2);
   pars(196) = me2(0,0);
   pars(197) = me2(0,1);
   pars(198) = me2(0,2);
   pars(199) = me2(1,0);
   pars(200) = me2(1,1);
   pars(201) = me2(1,2);
   pars(202) = me2(2,0);
   pars(203) = me2(2,1);
   pars(204) = me2(2,2);
   pars(205) = ms2;
   pars(206) = msbar2;
   pars(207) = mH1I2(0,0);
   pars(208) = mH1I2(0,1);
   pars(209) = mH1I2(1,0);
   pars(210) = mH1I2(1,1);
   pars(211) = mH2I2(0,0);
   pars(212) = mH2I2(0,1);
   pars(213) = mH2I2(1,0);
   pars(214) = mH2I2(1,1);
   pars(215) = mSI2(0,0);
   pars(216) = mSI2(0,1);
   pars(217) = mSI2(0,2);
   pars(218) = mSI2(1,0);
   pars(219) = mSI2(1,1);
   pars(220) = mSI2(1,2);
   pars(221) = mSI2(2,0);
   pars(222) = mSI2(2,1);
   pars(223) = mSI2(2,2);
   pars(224) = mDx2(0,0);
   pars(225) = mDx2(0,1);
   pars(226) = mDx2(0,2);
   pars(227) = mDx2(1,0);
   pars(228) = mDx2(1,1);
   pars(229) = mDx2(1,2);
   pars(230) = mDx2(2,0);
   pars(231) = mDx2(2,1);
   pars(232) = mDx2(2,2);
   pars(233) = mDxbar2(0,0);
   pars(234) = mDxbar2(0,1);
   pars(235) = mDxbar2(0,2);
   pars(236) = mDxbar2(1,0);
   pars(237) = mDxbar2(1,1);
   pars(238) = mDxbar2(1,2);
   pars(239) = mDxbar2(2,0);
   pars(240) = mDxbar2(2,1);
   pars(241) = mDxbar2(2,2);
   pars(242) = mHp2;
   pars(243) = mHpbar2;
   pars(244) = mphi2;
   pars(245) = MassB;
   pars(246) = MassWB;
   pars(247) = MassG;
   pars(248) = MassBp;


   return pars;
}

void CSE6SSM_soft_parameters::print(std::ostream& ostr) const
{
   CSE6SSM_susy_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "soft parameters:\n"
           "----------------------------------------\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "ThE = " << ThE << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TSigmaL = " << TSigmaL << '\n';
   ostr << "TKappaPr = " << TKappaPr << '\n';
   ostr << "TSigmax = " << TSigmax << '\n';
   ostr << "TgD = " << TgD << '\n';
   ostr << "TKappa = " << TKappa << '\n';
   ostr << "TLambda12 = " << TLambda12 << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "Tfu = " << Tfu << '\n';
   ostr << "Tfd = " << Tfd << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMuPr = " << BMuPr << '\n';
   ostr << "BMuPhi = " << BMuPhi << '\n';
   ostr << "LXiF = " << LXiF << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "msbar2 = " << msbar2 << '\n';
   ostr << "mH1I2 = " << mH1I2 << '\n';
   ostr << "mH2I2 = " << mH2I2 << '\n';
   ostr << "mSI2 = " << mSI2 << '\n';
   ostr << "mDx2 = " << mDx2 << '\n';
   ostr << "mDxbar2 = " << mDxbar2 << '\n';
   ostr << "mHp2 = " << mHp2 << '\n';
   ostr << "mHpbar2 = " << mHpbar2 << '\n';
   ostr << "mphi2 = " << mphi2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassBp = " << MassBp << '\n';

}

void CSE6SSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   CSE6SSM_susy_parameters::set(pars);

   TYd(0,0) = pars(84);
   TYd(0,1) = pars(85);
   TYd(0,2) = pars(86);
   TYd(1,0) = pars(87);
   TYd(1,1) = pars(88);
   TYd(1,2) = pars(89);
   TYd(2,0) = pars(90);
   TYd(2,1) = pars(91);
   TYd(2,2) = pars(92);
   ThE(0,0) = pars(93);
   ThE(0,1) = pars(94);
   ThE(1,0) = pars(95);
   ThE(1,1) = pars(96);
   ThE(2,0) = pars(97);
   ThE(2,1) = pars(98);
   TYe(0,0) = pars(99);
   TYe(0,1) = pars(100);
   TYe(0,2) = pars(101);
   TYe(1,0) = pars(102);
   TYe(1,1) = pars(103);
   TYe(1,2) = pars(104);
   TYe(2,0) = pars(105);
   TYe(2,1) = pars(106);
   TYe(2,2) = pars(107);
   TSigmaL = pars(108);
   TKappaPr = pars(109);
   TSigmax = pars(110);
   TgD(0,0) = pars(111);
   TgD(0,1) = pars(112);
   TgD(0,2) = pars(113);
   TgD(1,0) = pars(114);
   TgD(1,1) = pars(115);
   TgD(1,2) = pars(116);
   TgD(2,0) = pars(117);
   TgD(2,1) = pars(118);
   TgD(2,2) = pars(119);
   TKappa(0,0) = pars(120);
   TKappa(0,1) = pars(121);
   TKappa(0,2) = pars(122);
   TKappa(1,0) = pars(123);
   TKappa(1,1) = pars(124);
   TKappa(1,2) = pars(125);
   TKappa(2,0) = pars(126);
   TKappa(2,1) = pars(127);
   TKappa(2,2) = pars(128);
   TLambda12(0,0) = pars(129);
   TLambda12(0,1) = pars(130);
   TLambda12(1,0) = pars(131);
   TLambda12(1,1) = pars(132);
   TLambdax = pars(133);
   Tfu(0,0) = pars(134);
   Tfu(0,1) = pars(135);
   Tfu(1,0) = pars(136);
   Tfu(1,1) = pars(137);
   Tfu(2,0) = pars(138);
   Tfu(2,1) = pars(139);
   Tfd(0,0) = pars(140);
   Tfd(0,1) = pars(141);
   Tfd(1,0) = pars(142);
   Tfd(1,1) = pars(143);
   Tfd(2,0) = pars(144);
   Tfd(2,1) = pars(145);
   TYu(0,0) = pars(146);
   TYu(0,1) = pars(147);
   TYu(0,2) = pars(148);
   TYu(1,0) = pars(149);
   TYu(1,1) = pars(150);
   TYu(1,2) = pars(151);
   TYu(2,0) = pars(152);
   TYu(2,1) = pars(153);
   TYu(2,2) = pars(154);
   BMuPr = pars(155);
   BMuPhi = pars(156);
   LXiF = pars(157);
   mq2(0,0) = pars(158);
   mq2(0,1) = pars(159);
   mq2(0,2) = pars(160);
   mq2(1,0) = pars(161);
   mq2(1,1) = pars(162);
   mq2(1,2) = pars(163);
   mq2(2,0) = pars(164);
   mq2(2,1) = pars(165);
   mq2(2,2) = pars(166);
   ml2(0,0) = pars(167);
   ml2(0,1) = pars(168);
   ml2(0,2) = pars(169);
   ml2(1,0) = pars(170);
   ml2(1,1) = pars(171);
   ml2(1,2) = pars(172);
   ml2(2,0) = pars(173);
   ml2(2,1) = pars(174);
   ml2(2,2) = pars(175);
   mHd2 = pars(176);
   mHu2 = pars(177);
   md2(0,0) = pars(178);
   md2(0,1) = pars(179);
   md2(0,2) = pars(180);
   md2(1,0) = pars(181);
   md2(1,1) = pars(182);
   md2(1,2) = pars(183);
   md2(2,0) = pars(184);
   md2(2,1) = pars(185);
   md2(2,2) = pars(186);
   mu2(0,0) = pars(187);
   mu2(0,1) = pars(188);
   mu2(0,2) = pars(189);
   mu2(1,0) = pars(190);
   mu2(1,1) = pars(191);
   mu2(1,2) = pars(192);
   mu2(2,0) = pars(193);
   mu2(2,1) = pars(194);
   mu2(2,2) = pars(195);
   me2(0,0) = pars(196);
   me2(0,1) = pars(197);
   me2(0,2) = pars(198);
   me2(1,0) = pars(199);
   me2(1,1) = pars(200);
   me2(1,2) = pars(201);
   me2(2,0) = pars(202);
   me2(2,1) = pars(203);
   me2(2,2) = pars(204);
   ms2 = pars(205);
   msbar2 = pars(206);
   mH1I2(0,0) = pars(207);
   mH1I2(0,1) = pars(208);
   mH1I2(1,0) = pars(209);
   mH1I2(1,1) = pars(210);
   mH2I2(0,0) = pars(211);
   mH2I2(0,1) = pars(212);
   mH2I2(1,0) = pars(213);
   mH2I2(1,1) = pars(214);
   mSI2(0,0) = pars(215);
   mSI2(0,1) = pars(216);
   mSI2(0,2) = pars(217);
   mSI2(1,0) = pars(218);
   mSI2(1,1) = pars(219);
   mSI2(1,2) = pars(220);
   mSI2(2,0) = pars(221);
   mSI2(2,1) = pars(222);
   mSI2(2,2) = pars(223);
   mDx2(0,0) = pars(224);
   mDx2(0,1) = pars(225);
   mDx2(0,2) = pars(226);
   mDx2(1,0) = pars(227);
   mDx2(1,1) = pars(228);
   mDx2(1,2) = pars(229);
   mDx2(2,0) = pars(230);
   mDx2(2,1) = pars(231);
   mDx2(2,2) = pars(232);
   mDxbar2(0,0) = pars(233);
   mDxbar2(0,1) = pars(234);
   mDxbar2(0,2) = pars(235);
   mDxbar2(1,0) = pars(236);
   mDxbar2(1,1) = pars(237);
   mDxbar2(1,2) = pars(238);
   mDxbar2(2,0) = pars(239);
   mDxbar2(2,1) = pars(240);
   mDxbar2(2,2) = pars(241);
   mHp2 = pars(242);
   mHpbar2 = pars(243);
   mphi2 = pars(244);
   MassB = pars(245);
   MassWB = pars(246);
   MassG = pars(247);
   MassBp = pars(248);

}

void CSE6SSM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   TRACE_STRUCT.traceAdjfdTfd = Re((fd.adjoint()*Tfd).trace());
   TRACE_STRUCT.traceAdjYdTYd = Re((Yd.adjoint()*TYd).trace());
   TRACE_STRUCT.traceAdjYeTYe = Re((Ye.adjoint()*TYe).trace());
   TRACE_STRUCT.tracefdAdjfd = Re((fd*fd.adjoint()).trace());
   TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
   TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfu = Re((fu*fu.adjoint()).trace());
   TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappa = Re((Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12 = Re((Lambda12*(Lambda12).adjoint())
      .trace());
   TRACE_STRUCT.traceAdjfuTfu = Re((fu.adjoint()*Tfu).trace());
   TRACE_STRUCT.traceAdjYuTYu = Re((Yu.adjoint()*TYu).trace());
   TRACE_STRUCT.traceAdjKappaTKappa = Re(((Kappa).adjoint()*TKappa).trace());
   TRACE_STRUCT.traceAdjLambda12TLambda12 = Re(((Lambda12).adjoint()*TLambda12)
      .trace());
   TRACE_STRUCT.tracefdAdjfdTfdAdjfd = Re((fd*fd.adjoint()*Tfd*fd.adjoint())
      .trace());
   TRACE_STRUCT.tracefdAdjfdTfuAdjfu = Re((fd*fd.adjoint()*Tfu*fu.adjoint())
      .trace());
   TRACE_STRUCT.tracefuAdjfuTfdAdjfd = Re((fu*fu.adjoint()*Tfd*fd.adjoint())
      .trace());
   TRACE_STRUCT.tracehEAdjhETYeAdjYe = Re((hE*hE.adjoint()*TYe*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYdTYdAdjYd = Re((Yd*Yd.adjoint()*TYd*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYuTYuAdjYd = Re((Yd*Yu.adjoint()*TYu*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYeThEAdjhE = Re((Ye*Ye.adjoint()*ThE*hE.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYeTYeAdjYe = Re((Ye*Ye.adjoint()*TYe*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYdTYdAdjYu = Re((Yu*Yd.adjoint()*TYd*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceAdjfdTfdconjLambda12TpLambda12 = Re((fd.adjoint()*Tfd*
      Lambda12.conjugate()*(Lambda12).transpose()).trace());
   TRACE_STRUCT.traceAdjgDTpYdconjYdTgD = Re((gD.adjoint()*Yd.transpose()*
      Yd.conjugate()*TgD).trace());
   TRACE_STRUCT.traceAdjYdTYdconjgDTpgD = Re((Yd.adjoint()*TYd*gD.conjugate()*
      gD.transpose()).trace());
   TRACE_STRUCT.traceAdjLambda12TpfdconjfdTLambda12 = Re(((Lambda12).adjoint()*
      fd.transpose()*fd.conjugate()*TLambda12).trace());
   TRACE_STRUCT.traceAdjgDTgD = Re((gD.adjoint()*TgD).trace());
   TRACE_STRUCT.traceAdjhEThE = Re((hE.adjoint()*ThE).trace());
   TRACE_STRUCT.tracegDAdjgD = Re((gD*gD.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhE = Re((hE*hE.adjoint()).trace());
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
   TRACE_STRUCT.tracefuAdjhEThEAdjfu = Re((fu*hE.adjoint()*ThE*fu.adjoint())
      .trace());
   TRACE_STRUCT.tracegDAdjgDTgDAdjgD = Re((gD*gD.adjoint()*TgD*gD.adjoint())
      .trace());
   TRACE_STRUCT.tracegDAdjKappaTKappaAdjgD = Re((gD*(Kappa).adjoint()*TKappa*
      gD.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjfuTfuAdjhE = Re((hE*fu.adjoint()*Tfu*hE.adjoint())
      .trace());
   TRACE_STRUCT.tracehEAdjhEThEAdjhE = Re((hE*hE.adjoint()*ThE*hE.adjoint())
      .trace());
   TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjhE = Re((hE*(Lambda12).adjoint()*
      TLambda12*hE.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjgDTgDAdjKappa = Re((Kappa*gD.adjoint()*TgD*(Kappa)
      .adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjhEThEAdjLambda12 = Re((Lambda12*hE.adjoint()*
      ThE*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.traceAdjgDTpYuconjYuTgD = Re((gD.adjoint()*Yu.transpose()*
      Yu.conjugate()*TgD).trace());
   TRACE_STRUCT.traceAdjYuTYuconjgDTpgD = Re((Yu.adjoint()*TYu*gD.conjugate()*
      gD.transpose()).trace());
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
   TRACE_STRUCT.tracefuAdjLambda12TLambda12Adjfu = Re((fu*(Lambda12).adjoint()*
      TLambda12*fu.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = Re((Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa = Re((Kappa*(Kappa).adjoint()*
      TKappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjfuTfuAdjLambda12 = Re((Lambda12*fu.adjoint()*
      Tfu*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12 = Re((Lambda12*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12 = Re((Lambda12*(
      Lambda12).adjoint()*TLambda12*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfufuAdjfu = Re((fu*fu.adjoint()*fu*fu.adjoint())
      .trace());
   TRACE_STRUCT.tracefuAdjfuTfuAdjfu = Re((fu*fu.adjoint()*Tfu*fu.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYuTYuAdjYu = Re((Yu*Yu.adjoint()*TYu*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceconjTgDTpgD = Re((TgD.conjugate()*gD.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTgDTpTgD = Re((TgD.conjugate()*(TgD).transpose())
      .trace());
   TRACE_STRUCT.traceconjThETphE = Re((ThE.conjugate()*hE.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjThETpThE = Re((ThE.conjugate()*(ThE).transpose())
      .trace());
   TRACE_STRUCT.tracegDmDxbar2AdjgD = Re((gD*mDxbar2*gD.adjoint()).trace());
   TRACE_STRUCT.tracegDAdjgDmq2 = Re((gD*gD.adjoint()*mq2).trace());
   TRACE_STRUCT.tracehEAdjhEconjme2 = Re((hE*hE.adjoint()*me2.conjugate())
      .trace());
   TRACE_STRUCT.tracehEconjmH1I2AdjhE = Re((hE*mH1I2.conjugate()*hE.adjoint())
      .trace());
   TRACE_STRUCT.traceconjTfdTpTfd = Re((Tfd.conjugate()*(Tfd).transpose())
      .trace());
   TRACE_STRUCT.traceconjTYdTpTYd = Re((TYd.conjugate()*(TYd).transpose())
      .trace());
   TRACE_STRUCT.traceconjTYeTpTYe = Re((TYe.conjugate()*(TYe).transpose())
      .trace());
   TRACE_STRUCT.tracefdmH2I2Adjfd = Re((fd*mH2I2*fd.adjoint()).trace());
   TRACE_STRUCT.tracefdAdjfdconjmSI2 = Re((fd*fd.adjoint()*mSI2.conjugate())
      .trace());
   TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
   TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
   TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
   TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
   TRACE_STRUCT.traceconjTfdTpfd = Re((Tfd.conjugate()*fd.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYdTpYd = Re((TYd.conjugate()*Yd.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYeTpYe = Re((TYe.conjugate()*Ye.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTfuTpTfu = Re((Tfu.conjugate()*(Tfu).transpose())
      .trace());
   TRACE_STRUCT.traceconjTYuTpTYu = Re((TYu.conjugate()*(TYu).transpose())
      .trace());
   TRACE_STRUCT.tracefumH1I2Adjfu = Re((fu*mH1I2*fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfuconjmSI2 = Re((fu*fu.adjoint()*mSI2.conjugate())
      .trace());
   TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
   TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());
   TRACE_STRUCT.traceconjTfuTpfu = Re((Tfu.conjugate()*fu.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYuTpYu = Re((TYu.conjugate()*Yu.transpose()).trace())
      ;
   TRACE_STRUCT.tracegDAdjgDconjmq2 = Re((gD*gD.adjoint()*mq2.conjugate())
      .trace());
   TRACE_STRUCT.tracegDconjmDxbar2AdjgD = Re((gD*mDxbar2.conjugate()*gD.adjoint
      ()).trace());
   TRACE_STRUCT.tracehEmH1I2AdjhE = Re((hE*mH1I2*hE.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhEme2 = Re((hE*hE.adjoint()*me2).trace());
   TRACE_STRUCT.traceconjTKappaTpKappa = Re((TKappa.conjugate()*(Kappa)
      .transpose()).trace());
   TRACE_STRUCT.traceconjTKappaTpTKappa = Re((TKappa.conjugate()*(TKappa)
      .transpose()).trace());
   TRACE_STRUCT.traceconjTLambda12TpLambda12 = Re((TLambda12.conjugate()*(
      Lambda12).transpose()).trace());
   TRACE_STRUCT.traceconjTLambda12TpTLambda12 = Re((TLambda12.conjugate()*(
      TLambda12).transpose()).trace());
   TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12 = Re((mH1I2*(Lambda12).adjoint()*
      Lambda12).trace());
   TRACE_STRUCT.traceKappaAdjKappaconjmDx2 = Re((Kappa*(Kappa).adjoint()*
      mDx2.conjugate()).trace());
   TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa = Re((Kappa*mDxbar2.conjugate()*(
      Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2 = Re((Lambda12*(Lambda12)
      .adjoint()*mH2I2.conjugate()).trace());
   TRACE_STRUCT.tracefdAdjfdTfdAdjTfd = Re((fd*fd.adjoint()*Tfd*(Tfd).adjoint()
      ).trace());
   TRACE_STRUCT.tracefdAdjfdTfuAdjTfu = Re((fd*fd.adjoint()*Tfu*(Tfu).adjoint()
      ).trace());
   TRACE_STRUCT.tracefdAdjTfdTfdAdjfd = Re((fd*(Tfd).adjoint()*Tfd*fd.adjoint()
      ).trace());
   TRACE_STRUCT.tracefdAdjTfdTfuAdjfu = Re((fd*(Tfd).adjoint()*Tfu*fu.adjoint()
      ).trace());
   TRACE_STRUCT.tracefdconjTLambda12TpTLambda12Adjfd = Re((fd*
      TLambda12.conjugate()*(TLambda12).transpose()*fd.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfuTfdAdjTfd = Re((fu*fu.adjoint()*Tfd*(Tfd).adjoint()
      ).trace());
   TRACE_STRUCT.tracefuAdjTfuTfdAdjfd = Re((fu*(Tfu).adjoint()*Tfd*fd.adjoint()
      ).trace());
   TRACE_STRUCT.tracegDAdjgDTpTYdconjTYd = Re((gD*gD.adjoint()*(TYd).transpose(
      )*TYd.conjugate()).trace());
   TRACE_STRUCT.tracehEAdjhETYeAdjTYe = Re((hE*hE.adjoint()*TYe*(TYe).adjoint()
      ).trace());
   TRACE_STRUCT.tracehEAdjThETYeAdjYe = Re((hE*(ThE).adjoint()*TYe*Ye.adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = Re((Yd*Yd.adjoint()*TYd*(TYd).adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = Re((Yd*Yu.adjoint()*TYu*(TYd).adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = Re((Yd*(TYd).adjoint()*TYd*Yd.adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = Re((Yd*(TYu).adjoint()*TYu*Yd.adjoint()
      ).trace());
   TRACE_STRUCT.traceYdconjTgDTpTgDAdjYd = Re((Yd*TgD.conjugate()*(TgD)
      .transpose()*Yd.adjoint()).trace());
   TRACE_STRUCT.traceYeAdjYeThEAdjThE = Re((Ye*Ye.adjoint()*ThE*(ThE).adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = Re((Ye*Ye.adjoint()*TYe*(TYe).adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjTYeThEAdjhE = Re((Ye*(TYe).adjoint()*ThE*hE.adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = Re((Ye*(TYe).adjoint()*TYe*Ye.adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = Re((Yu*Yd.adjoint()*TYd*(TYu).adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = Re((Yu*(TYd).adjoint()*TYd*Yu.adjoint()
      ).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12TpTfdconjTfd = Re((Lambda12*(Lambda12)
      .adjoint()*(Tfd).transpose()*Tfd.conjugate()).trace());
   TRACE_STRUCT.traceAdjfdTfdconjTLambda12TpLambda12 = Re((fd.adjoint()*Tfd*
      TLambda12.conjugate()*(Lambda12).transpose()).trace());
   TRACE_STRUCT.traceAdjgDTpYdconjTYdTgD = Re((gD.adjoint()*Yd.transpose()*
      TYd.conjugate()*TgD).trace());
   TRACE_STRUCT.traceAdjYdTYdconjTgDTpgD = Re((Yd.adjoint()*TYd*TgD.conjugate()
      *gD.transpose()).trace());
   TRACE_STRUCT.traceAdjLambda12TpfdconjTfdTLambda12 = Re(((Lambda12).adjoint()
      *fd.transpose()*Tfd.conjugate()*TLambda12).trace());
   TRACE_STRUCT.tracefdmH2I2AdjfdfdAdjfd = Re((fd*mH2I2*fd.adjoint()*fd*
      fd.adjoint()).trace());
   TRACE_STRUCT.tracefdmH2I2AdjfdfuAdjfu = Re((fd*mH2I2*fd.adjoint()*fu*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefdAdjfdfdmH2I2Adjfd = Re((fd*fd.adjoint()*fd*mH2I2*
      fd.adjoint()).trace());
   TRACE_STRUCT.tracefdAdjfdfumH1I2Adjfu = Re((fd*fd.adjoint()*fu*mH1I2*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefdAdjfdfuAdjfuconjmSI2 = Re((fd*fd.adjoint()*fu*fu.adjoint
      ()*mSI2.conjugate()).trace());
   TRACE_STRUCT.tracefdAdjfdconjmSI2fdAdjfd = Re((fd*fd.adjoint()*
      mSI2.conjugate()*fd*fd.adjoint()).trace());
   TRACE_STRUCT.tracefdAdjfdconjmSI2fuAdjfu = Re((fd*fd.adjoint()*
      mSI2.conjugate()*fu*fu.adjoint()).trace());
   TRACE_STRUCT.tracefdconjLambda12TpLambda12AdjfdconjmSI2 = Re((fd*
      Lambda12.conjugate()*(Lambda12).transpose()*fd.adjoint()*mSI2.conjugate())
      .trace());
   TRACE_STRUCT.tracegDAdjgDconjmq2TpYdconjYd = Re((gD*gD.adjoint()*
      mq2.conjugate()*Yd.transpose()*Yd.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjgDTpYdconjmd2conjYd = Re((gD*gD.adjoint()*
      Yd.transpose()*md2.conjugate()*Yd.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjgDTpYdconjYdconjmq2 = Re((gD*gD.adjoint()*
      Yd.transpose()*Yd.conjugate()*mq2.conjugate()).trace());
   TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYdconjYd = Re((gD*mDxbar2.conjugate()*
      gD.adjoint()*Yd.transpose()*Yd.conjugate()).trace());
   TRACE_STRUCT.tracehEmH1I2AdjhEYeAdjYe = Re((hE*mH1I2*hE.adjoint()*Ye*
      Ye.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhEme2YeAdjYe = Re((hE*hE.adjoint()*me2*Ye*Ye.adjoint(
      )).trace());
   TRACE_STRUCT.tracehEAdjhEYeml2AdjYe = Re((hE*hE.adjoint()*Ye*ml2*Ye.adjoint(
      )).trace());
   TRACE_STRUCT.tracehEAdjhEYeAdjYeme2 = Re((hE*hE.adjoint()*Ye*Ye.adjoint()*
      me2).trace());
   TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = Re((md2*Yd*Yd.adjoint()*Yd*Yd.adjoint(
      )).trace());
   TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = Re((md2*Yd*Yu.adjoint()*Yu*Yd.adjoint(
      )).trace());
   TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = Re((me2*Ye*Ye.adjoint()*Ye*Ye.adjoint(
      )).trace());
   TRACE_STRUCT.tracemH1I2AdjLambda12TpfdconjfdLambda12 = Re((mH1I2*(Lambda12)
      .adjoint()*fd.transpose()*fd.conjugate()*Lambda12).trace());
   TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = Re((ml2*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye).trace());
   TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = Re((mq2*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd).trace());
   TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = Re((mq2*Yd.adjoint()*Yd*Yu.adjoint()*
      Yu).trace());
   TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = Re((mq2*Yu.adjoint()*Yu*Yd.adjoint()*
      Yd).trace());
   TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = Re((mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint(
      )).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Tpfdconjfd = Re((Lambda12*(
      Lambda12).adjoint()*mH2I2.conjugate()*fd.transpose()*fd.conjugate()).trace()
      );
   TRACE_STRUCT.traceLambda12AdjLambda12TpfdconjfdconjmH2I2 = Re((Lambda12*(
      Lambda12).adjoint()*fd.transpose()*fd.conjugate()*mH2I2.conjugate()).trace()
      );
   TRACE_STRUCT.tracefuAdjfuTfuAdjTfu = Re((fu*fu.adjoint()*Tfu*(Tfu).adjoint()
      ).trace());
   TRACE_STRUCT.tracefuAdjhEThEAdjTfu = Re((fu*hE.adjoint()*ThE*(Tfu).adjoint()
      ).trace());
   TRACE_STRUCT.tracefuAdjLambda12TLambda12AdjTfu = Re((fu*(Lambda12).adjoint()
      *TLambda12*(Tfu).adjoint()).trace());
   TRACE_STRUCT.tracefuAdjTfuTfuAdjfu = Re((fu*(Tfu).adjoint()*Tfu*fu.adjoint()
      ).trace());
   TRACE_STRUCT.tracefuAdjThEThEAdjfu = Re((fu*(ThE).adjoint()*ThE*fu.adjoint()
      ).trace());
   TRACE_STRUCT.tracefuAdjTLambda12TLambda12Adjfu = Re((fu*(TLambda12).adjoint(
      )*TLambda12*fu.adjoint()).trace());
   TRACE_STRUCT.tracegDAdjgDTpTYuconjTYu = Re((gD*gD.adjoint()*(TYu).transpose(
      )*TYu.conjugate()).trace());
   TRACE_STRUCT.tracehEAdjfuTfuAdjThE = Re((hE*fu.adjoint()*Tfu*(ThE).adjoint()
      ).trace());
   TRACE_STRUCT.tracehEAdjTfuTfuAdjhE = Re((hE*(Tfu).adjoint()*Tfu*hE.adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = Re((Yu*Yu.adjoint()*TYu*(TYu).adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = Re((Yu*(TYu).adjoint()*TYu*Yu.adjoint()
      ).trace());
   TRACE_STRUCT.traceYuconjTgDTpTgDAdjYu = Re((Yu*TgD.conjugate()*(TgD)
      .transpose()*Yu.adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjfuTfuAdjTLambda12 = Re((Lambda12*fu.adjoint()*
      Tfu*(TLambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjTfuTfuAdjLambda12 = Re((Lambda12*(Tfu).adjoint(
      )*Tfu*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.traceAdjgDTpYuconjTYuTgD = Re((gD.adjoint()*Yu.transpose()*
      TYu.conjugate()*TgD).trace());
   TRACE_STRUCT.traceAdjYuTYuconjTgDTpgD = Re((Yu.adjoint()*TYu*TgD.conjugate()
      *gD.transpose()).trace());
   TRACE_STRUCT.tracefumH1I2AdjfufuAdjfu = Re((fu*mH1I2*fu.adjoint()*fu*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefumH1I2AdjhEhEAdjfu = Re((fu*mH1I2*hE.adjoint()*hE*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefumH1I2AdjLambda12Lambda12Adjfu = Re((fu*mH1I2*(Lambda12)
      .adjoint()*Lambda12*fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfufumH1I2Adjfu = Re((fu*fu.adjoint()*fu*mH1I2*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjfuconjmSI2fuAdjfu = Re((fu*fu.adjoint()*
      mSI2.conjugate()*fu*fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjhEhEmH1I2Adjfu = Re((fu*hE.adjoint()*hE*mH1I2*
      fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjhEhEAdjfuconjmSI2 = Re((fu*hE.adjoint()*hE*fu.adjoint
      ()*mSI2.conjugate()).trace());
   TRACE_STRUCT.tracefuAdjhEme2hEAdjfu = Re((fu*hE.adjoint()*me2*hE*fu.adjoint(
      )).trace());
   TRACE_STRUCT.tracefuAdjLambda12Lambda12mH1I2Adjfu = Re((fu*(Lambda12)
      .adjoint()*Lambda12*mH1I2*fu.adjoint()).trace());
   TRACE_STRUCT.tracefuAdjLambda12Lambda12AdjfuconjmSI2 = Re((fu*(Lambda12)
      .adjoint()*Lambda12*fu.adjoint()*mSI2.conjugate()).trace());
   TRACE_STRUCT.tracefuAdjLambda12conjmH2I2Lambda12Adjfu = Re((fu*(Lambda12)
      .adjoint()*mH2I2.conjugate()*Lambda12*fu.adjoint()).trace());
   TRACE_STRUCT.tracegDAdjgDconjmq2TpYuconjYu = Re((gD*gD.adjoint()*
      mq2.conjugate()*Yu.transpose()*Yu.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjgDTpYuconjmu2conjYu = Re((gD*gD.adjoint()*
      Yu.transpose()*mu2.conjugate()*Yu.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjgDTpYuconjYuconjmq2 = Re((gD*gD.adjoint()*
      Yu.transpose()*Yu.conjugate()*mq2.conjugate()).trace());
   TRACE_STRUCT.tracegDconjmDxbar2AdjgDTpYuconjYu = Re((gD*mDxbar2.conjugate()*
      gD.adjoint()*Yu.transpose()*Yu.conjugate()).trace());
   TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = Re((mq2*Yu.adjoint()*Yu*Yu.adjoint()*
      Yu).trace());
   TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = Re((mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint(
      )).trace());
   TRACE_STRUCT.tracegDAdjKappaTKappaAdjTgD = Re((gD*(Kappa).adjoint()*TKappa*(
      TgD).adjoint()).trace());
   TRACE_STRUCT.tracegDAdjTKappaTKappaAdjgD = Re((gD*(TKappa).adjoint()*TKappa*
      gD.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjLambda12TLambda12AdjThE = Re((hE*(Lambda12).adjoint()
      *TLambda12*(ThE).adjoint()).trace());
   TRACE_STRUCT.tracehEAdjTLambda12TLambda12AdjhE = Re((hE*(TLambda12).adjoint(
      )*TLambda12*hE.adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjgDTgDAdjTKappa = Re((Kappa*gD.adjoint()*TgD*(
      TKappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa = Re((Kappa*(Kappa).adjoint()
      *TKappa*(TKappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjTgDTgDAdjKappa = Re((Kappa*(TgD).adjoint()*TgD*(
      Kappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa = Re((Kappa*(TKappa).adjoint(
      )*TKappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjhEThEAdjTLambda12 = Re((Lambda12*hE.adjoint()*
      ThE*(TLambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12 = Re((Lambda12*(
      Lambda12).adjoint()*TLambda12*(TLambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjThEThEAdjLambda12 = Re((Lambda12*(ThE).adjoint(
      )*ThE*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12 = Re((Lambda12*(
      TLambda12).adjoint()*TLambda12*(Lambda12).adjoint()).trace());
   TRACE_STRUCT.tracegDAdjKappaKappaAdjgDconjmq2 = Re((gD*(Kappa).adjoint()*
      Kappa*gD.adjoint()*mq2.conjugate()).trace());
   TRACE_STRUCT.tracegDAdjKappaKappaconjmDxbar2AdjgD = Re((gD*(Kappa).adjoint()
      *Kappa*mDxbar2.conjugate()*gD.adjoint()).trace());
   TRACE_STRUCT.tracegDAdjKappaconjmDx2KappaAdjgD = Re((gD*(Kappa).adjoint()*
      mDx2.conjugate()*Kappa*gD.adjoint()).trace());
   TRACE_STRUCT.tracegDconjmDxbar2AdjKappaKappaAdjgD = Re((gD*mDxbar2.conjugate
      ()*(Kappa).adjoint()*Kappa*gD.adjoint()).trace());
   TRACE_STRUCT.tracehEmH1I2AdjLambda12Lambda12AdjhE = Re((hE*mH1I2*(Lambda12)
      .adjoint()*Lambda12*hE.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjLambda12Lambda12mH1I2AdjhE = Re((hE*(Lambda12)
      .adjoint()*Lambda12*mH1I2*hE.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjLambda12Lambda12AdjhEme2 = Re((hE*(Lambda12).adjoint(
      )*Lambda12*hE.adjoint()*me2).trace());
   TRACE_STRUCT.tracehEAdjLambda12conjmH2I2Lambda12AdjhE = Re((hE*(Lambda12)
      .adjoint()*mH2I2.conjugate()*Lambda12*hE.adjoint()).trace());
   TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 = Re((mH1I2*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12).trace());
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2 = Re((Kappa*(Kappa)
      .adjoint()*Kappa*(Kappa).adjoint()*mDx2.conjugate()).trace());
   TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa = Re((Kappa*(Kappa)
      .adjoint()*Kappa*mDxbar2.conjugate()*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa = Re((Kappa*(Kappa)
      .adjoint()*mDx2.conjugate()*Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa = Re((Kappa*
      mDxbar2.conjugate()*(Kappa).adjoint()*Kappa*(Kappa).adjoint()).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 = Re((
      Lambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*mH2I2.conjugate(
      )).trace());
   TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 = Re((
      Lambda12*(Lambda12).adjoint()*mH2I2.conjugate()*Lambda12*(Lambda12).adjoint(
      )).trace());
   TRACE_STRUCT.tracegDAdjgDTgDAdjTgD = Re((gD*gD.adjoint()*TgD*(TgD).adjoint()
      ).trace());
   TRACE_STRUCT.tracegDAdjTgDTgDAdjgD = Re((gD*(TgD).adjoint()*TgD*gD.adjoint()
      ).trace());
   TRACE_STRUCT.tracehEAdjhEThEAdjThE = Re((hE*hE.adjoint()*ThE*(ThE).adjoint()
      ).trace());
   TRACE_STRUCT.tracehEAdjThEThEAdjhE = Re((hE*(ThE).adjoint()*ThE*hE.adjoint()
      ).trace());
   TRACE_STRUCT.tracegDAdjgDgDconjmDxbar2AdjgD = Re((gD*gD.adjoint()*gD*
      mDxbar2.conjugate()*gD.adjoint()).trace());
   TRACE_STRUCT.tracegDAdjgDconjmq2gDAdjgD = Re((gD*gD.adjoint()*mq2.conjugate(
      )*gD*gD.adjoint()).trace());
   TRACE_STRUCT.tracegDconjmDxbar2AdjgDgDAdjgD = Re((gD*mDxbar2.conjugate()*
      gD.adjoint()*gD*gD.adjoint()).trace());
   TRACE_STRUCT.tracehEmH1I2AdjhEhEAdjhE = Re((hE*mH1I2*hE.adjoint()*hE*
      hE.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhEhEmH1I2AdjhE = Re((hE*hE.adjoint()*hE*mH1I2*
      hE.adjoint()).trace());
   TRACE_STRUCT.tracehEAdjhEhEAdjhEme2 = Re((hE*hE.adjoint()*hE*hE.adjoint()*
      me2).trace());
   TRACE_STRUCT.tracehEAdjhEme2hEAdjhE = Re((hE*hE.adjoint()*me2*hE*hE.adjoint(
      )).trace());
   TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 - mHp2 + mHpbar2 + mHu2
      + (md2).trace() - (mDx2).trace() + (mDxbar2).trace() + (me2).trace() - (
      mH1I2).trace() + (mH2I2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2)
      .trace()));
   TRACE_STRUCT.Tr14 = Re(0.15811388300841897*g1p*(-6*mHd2 + 4*mHp2 - 4*mHpbar2
      - 4*mHu2 + ms2*QS - msbar2*QS + 6*(md2).trace() - 6*(mDx2).trace() - 9*(
      mDxbar2).trace() + (me2).trace() - 6*(mH1I2).trace() - 4*(mH2I2).trace() + 4
      *(ml2).trace() + 6*(mq2).trace() + 5*(mSI2).trace() + 3*(mu2).trace()));
   TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHp2 + 3*mHpbar2 + 3*mHu2
      + 2*(md2).trace() + 2*(mDx2).trace() + 2*(mDxbar2).trace() + 6*(me2).trace()
      + 3*(mH1I2).trace() + 3*(mH2I2).trace() + 3*(ml2).trace() + (mq2).trace() +
      8*(mu2).trace()));
   TRACE_STRUCT.Tr2U114 = Re(0.1224744871391589*g1*g1p*(3*mHd2 - 2*mHp2 - 2*
      mHpbar2 - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace()
      + (me2).trace() + 3*(mH1I2).trace() - 2*(mH2I2).trace() - 2*(ml2).trace() +
      (mq2).trace() - 2*(mu2).trace()));
   TRACE_STRUCT.Tr31 = Re(0.006454972243679028*g1*(60*(mHd2 - mHu2)*AbsSqr(
      Lambdax) + 60*(mHp2 - mHpbar2)*AbsSqr(SigmaL) - 18*mHd2*Sqr(g1) - 18*mHp2*
      Sqr(g1) + 18*mHpbar2*Sqr(g1) + 18*mHu2*Sqr(g1) - 27*mHd2*Sqr(g1p) - 12*mHp2*
      Sqr(g1p) + 12*mHpbar2*Sqr(g1p) + 12*mHu2*Sqr(g1p) - 90*mHd2*Sqr(g2) - 90*
      mHp2*Sqr(g2) + 90*mHpbar2*Sqr(g2) + 90*mHu2*Sqr(g2) + 8*Sqr(g1)*(md2).trace(
      ) + 12*Sqr(g1p)*(md2).trace() + 160*Sqr(g3)*(md2).trace() - 8*Sqr(g1)*(mDx2)
      .trace() - 12*Sqr(g1p)*(mDx2).trace() - 160*Sqr(g3)*(mDx2).trace() + 8*Sqr(
      g1)*(mDxbar2).trace() + 27*Sqr(g1p)*(mDxbar2).trace() + 160*Sqr(g3)*(mDxbar2
      ).trace() + 72*Sqr(g1)*(me2).trace() + 3*Sqr(g1p)*(me2).trace() - 18*Sqr(g1)
      *(mH1I2).trace() - 27*Sqr(g1p)*(mH1I2).trace() - 90*Sqr(g2)*(mH1I2).trace()
      + 18*Sqr(g1)*(mH2I2).trace() + 12*Sqr(g1p)*(mH2I2).trace() + 90*Sqr(g2)*(
      mH2I2).trace() - 18*Sqr(g1)*(ml2).trace() - 12*Sqr(g1p)*(ml2).trace() - 90*
      Sqr(g2)*(ml2).trace() + 2*Sqr(g1)*(mq2).trace() + 3*Sqr(g1p)*(mq2).trace() +
      90*Sqr(g2)*(mq2).trace() + 160*Sqr(g3)*(mq2).trace() - 64*Sqr(g1)*(mu2)
      .trace() - 6*Sqr(g1p)*(mu2).trace() - 320*Sqr(g3)*(mu2).trace() + 60*mHd2*(
      fd*fd.adjoint()).trace() - 60*mHu2*(fu*fu.adjoint()).trace() + 180*mHp2*(gD*
      gD.adjoint()).trace() + 60*mHp2*(hE*hE.adjoint()).trace() + 180*mHd2*(Yd*
      Yd.adjoint()).trace() + 60*mHd2*(Ye*Ye.adjoint()).trace() - 180*mHu2*(Yu*
      Yu.adjoint()).trace() - 60*(fd*mH2I2.conjugate()*fd.adjoint()).trace() + 60*
      (fu*mH1I2.conjugate()*fu.adjoint()).trace() - 120*(gD*mDxbar2*gD.adjoint())
      .trace() - 60*(gD*gD.adjoint()*mq2).trace() - 120*(hE*hE.adjoint()*
      me2.conjugate()).trace() + 60*(hE*mH1I2.conjugate()*hE.adjoint()).trace() +
      60*(mDx2*Kappa*(Kappa).adjoint()).trace() - 60*(mDxbar2*(Kappa).adjoint()*
      Kappa).trace() - 60*(mH2I2*Lambda12*(Lambda12).adjoint()).trace() - 120*(Yd*
      Yd.adjoint()*md2.conjugate()).trace() - 60*(Yd*mq2.conjugate()*Yd.adjoint())
      .trace() - 120*(Ye*Ye.adjoint()*me2.conjugate()).trace() + 60*(Ye*
      ml2.conjugate()*Ye.adjoint()).trace() + 240*(Yu*Yu.adjoint()*mu2.conjugate()
      ).trace() - 60*(Yu*mq2.conjugate()*Yu.adjoint()).trace() + 60*(Lambda12*
      mH1I2.conjugate()*(Lambda12).adjoint()).trace()));
   TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHp2 + mHpbar2 + mHu2 + (mH1I2).trace() +
      (mH2I2).trace() + (ml2).trace() + 3*(mq2).trace()));
   TRACE_STRUCT.Tr23 = Re(0.5*((md2).trace() + (mDx2).trace() + (mDxbar2).trace
      () + 2*(mq2).trace() + (mu2).trace()));
   TRACE_STRUCT.Tr2U141 = Re(0.1224744871391589*g1*g1p*(3*mHd2 - 2*mHp2 - 2*
      mHpbar2 - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace()
      + (me2).trace() + 3*(mH1I2).trace() - 2*(mH2I2).trace() - 2*(ml2).trace() +
      (mq2).trace() - 2*(mu2).trace()));
   TRACE_STRUCT.Tr2U144 = Re(0.025*Sqr(g1p)*(18*mHd2 + 8*mHp2 + 8*mHpbar2 + 8*
      mHu2 + ms2*Sqr(QS) + msbar2*Sqr(QS) + 12*(md2).trace() + 12*(mDx2).trace() +
      27*(mDxbar2).trace() + (me2).trace() + 18*(mH1I2).trace() + 8*(mH2I2).trace
      () + 8*(ml2).trace() + 6*(mq2).trace() + 25*(mSI2).trace() + 3*(mu2).trace()
      ));
   TRACE_STRUCT.Tr34 = Re(-0.003952847075210475*g1p*(-40*(3*mHd2 + 2*mHu2 - ms2
      *QS)*AbsSqr(Lambdax) + 20*(ms2 - msbar2)*QS*AbsSqr(Sigmax) + 80*mHp2*AbsSqr(
      SigmaL) - 80*mHpbar2*AbsSqr(SigmaL) + 36*mHd2*Sqr(g1) - 24*mHp2*Sqr(g1) + 24
      *mHpbar2*Sqr(g1) + 24*mHu2*Sqr(g1) + 54*mHd2*Sqr(g1p) - 16*mHp2*Sqr(g1p) +
      16*mHpbar2*Sqr(g1p) + 16*mHu2*Sqr(g1p) - ms2*Power(QS,3)*Sqr(g1p) + msbar2*
      Power(QS,3)*Sqr(g1p) + 180*mHd2*Sqr(g2) - 120*mHp2*Sqr(g2) + 120*mHpbar2*Sqr
      (g2) + 120*mHu2*Sqr(g2) - 16*Sqr(g1)*(md2).trace() - 24*Sqr(g1p)*(md2).trace
      () - 320*Sqr(g3)*(md2).trace() + 16*Sqr(g1)*(mDx2).trace() + 24*Sqr(g1p)*(
      mDx2).trace() + 320*Sqr(g3)*(mDx2).trace() + 24*Sqr(g1)*(mDxbar2).trace() +
      81*Sqr(g1p)*(mDxbar2).trace() + 480*Sqr(g3)*(mDxbar2).trace() - 24*Sqr(g1)*(
      me2).trace() - Sqr(g1p)*(me2).trace() + 36*Sqr(g1)*(mH1I2).trace() + 54*Sqr(
      g1p)*(mH1I2).trace() + 180*Sqr(g2)*(mH1I2).trace() + 24*Sqr(g1)*(mH2I2)
      .trace() + 16*Sqr(g1p)*(mH2I2).trace() + 120*Sqr(g2)*(mH2I2).trace() - 24*
      Sqr(g1)*(ml2).trace() - 16*Sqr(g1p)*(ml2).trace() - 120*Sqr(g2)*(ml2).trace(
      ) - 4*Sqr(g1)*(mq2).trace() - 6*Sqr(g1p)*(mq2).trace() - 180*Sqr(g2)*(mq2)
      .trace() - 320*Sqr(g3)*(mq2).trace() - 125*Sqr(g1p)*(mSI2).trace() - 32*Sqr(
      g1)*(mu2).trace() - 3*Sqr(g1p)*(mu2).trace() - 160*Sqr(g3)*(mu2).trace() -
      120*mHd2*(fd*fd.adjoint()).trace() - 80*mHu2*(fu*fu.adjoint()).trace() + 240
      *mHp2*(gD*gD.adjoint()).trace() + 80*mHp2*(hE*hE.adjoint()).trace() - 360*
      mHd2*(Yd*Yd.adjoint()).trace() - 120*mHd2*(Ye*Ye.adjoint()).trace() - 240*
      mHu2*(Yu*Yu.adjoint()).trace() + 60*ms2*QS*(Kappa*(Kappa).adjoint()).trace()
      + 40*ms2*QS*(Lambda12*(Lambda12).adjoint()).trace() + 200*(fd*fd.adjoint()*
      mSI2).trace() - 80*(fd*mH2I2.conjugate()*fd.adjoint()).trace() + 200*(fu*
      fu.adjoint()*mSI2).trace() - 120*(fu*mH1I2.conjugate()*fu.adjoint()).trace()
      - 360*(gD*mDxbar2*gD.adjoint()).trace() + 120*(gD*gD.adjoint()*mq2).trace()
      + 40*(hE*hE.adjoint()*me2.conjugate()).trace() - 120*(hE*mH1I2.conjugate()*
      hE.adjoint()).trace() - 120*(mDx2*Kappa*(Kappa).adjoint()).trace() - 180*(
      mDxbar2*(Kappa).adjoint()*Kappa).trace() - 80*(mH2I2*Lambda12*(Lambda12)
      .adjoint()).trace() + 240*(Yd*Yd.adjoint()*md2.conjugate()).trace() + 120*(
      Yd*mq2.conjugate()*Yd.adjoint()).trace() + 40*(Ye*Ye.adjoint()*me2.conjugate
      ()).trace() + 80*(Ye*ml2.conjugate()*Ye.adjoint()).trace() + 120*(Yu*
      Yu.adjoint()*mu2.conjugate()).trace() + 120*(Yu*mq2.conjugate()*Yu.adjoint()
      ).trace() - 120*(Lambda12*mH1I2.conjugate()*(Lambda12).adjoint()).trace()));

}

std::ostream& operator<<(std::ostream& ostr, const CSE6SSM_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
