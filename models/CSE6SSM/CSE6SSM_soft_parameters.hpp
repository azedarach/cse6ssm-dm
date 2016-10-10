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

#ifndef CSE6SSM_soft_parameters_H
#define CSE6SSM_soft_parameters_H

#include "rge.h"
#include "CSE6SSM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class CSE6SSM_soft_parameters : public CSE6SSM_susy_parameters {
public:
   explicit CSE6SSM_soft_parameters();
   CSE6SSM_soft_parameters(const CSE6SSM_susy_parameters& , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,2>&
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
);
   virtual ~CSE6SSM_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual const Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   CSE6SSM_soft_parameters calc_beta() const;
   virtual void clear();

   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) { TYd = TYd_; }
   void set_TYd(int i, int k, double value) { TYd(i,k) = value; }
   void set_ThE(const Eigen::Matrix<double,3,2>& ThE_) { ThE = ThE_; }
   void set_ThE(int i, int k, double value) { ThE(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) { TYe = TYe_; }
   void set_TYe(int i, int k, double value) { TYe(i,k) = value; }
   void set_TSigmaL(double TSigmaL_) { TSigmaL = TSigmaL_; }
   void set_TKappaPr(double TKappaPr_) { TKappaPr = TKappaPr_; }
   void set_TSigmax(double TSigmax_) { TSigmax = TSigmax_; }
   void set_TgD(const Eigen::Matrix<double,3,3>& TgD_) { TgD = TgD_; }
   void set_TgD(int i, int k, double value) { TgD(i,k) = value; }
   void set_TKappa(const Eigen::Matrix<double,3,3>& TKappa_) { TKappa = TKappa_; }
   void set_TKappa(int i, int k, double value) { TKappa(i,k) = value; }
   void set_TLambda12(const Eigen::Matrix<double,2,2>& TLambda12_) { TLambda12 = TLambda12_; }
   void set_TLambda12(int i, int k, double value) { TLambda12(i,k) = value; }
   void set_TLambdax(double TLambdax_) { TLambdax = TLambdax_; }
   void set_Tfu(const Eigen::Matrix<double,3,2>& Tfu_) { Tfu = Tfu_; }
   void set_Tfu(int i, int k, double value) { Tfu(i,k) = value; }
   void set_Tfd(const Eigen::Matrix<double,3,2>& Tfd_) { Tfd = Tfd_; }
   void set_Tfd(int i, int k, double value) { Tfd(i,k) = value; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) { TYu = TYu_; }
   void set_TYu(int i, int k, double value) { TYu(i,k) = value; }
   void set_BMuPr(double BMuPr_) { BMuPr = BMuPr_; }
   void set_BMuPhi(double BMuPhi_) { BMuPhi = BMuPhi_; }
   void set_LXiF(double LXiF_) { LXiF = LXiF_; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) { mq2 = mq2_; }
   void set_mq2(int i, int k, double value) { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) { ml2 = ml2_; }
   void set_ml2(int i, int k, double value) { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) { md2 = md2_; }
   void set_md2(int i, int k, double value) { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) { mu2 = mu2_; }
   void set_mu2(int i, int k, double value) { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) { me2 = me2_; }
   void set_me2(int i, int k, double value) { me2(i,k) = value; }
   void set_ms2(double ms2_) { ms2 = ms2_; }
   void set_msbar2(double msbar2_) { msbar2 = msbar2_; }
   void set_mH1I2(const Eigen::Matrix<double,2,2>& mH1I2_) { mH1I2 = mH1I2_; }
   void set_mH1I2(int i, int k, double value) { mH1I2(i,k) = value; }
   void set_mH2I2(const Eigen::Matrix<double,2,2>& mH2I2_) { mH2I2 = mH2I2_; }
   void set_mH2I2(int i, int k, double value) { mH2I2(i,k) = value; }
   void set_mSI2(const Eigen::Matrix<double,3,3>& mSI2_) { mSI2 = mSI2_; }
   void set_mSI2(int i, int k, double value) { mSI2(i,k) = value; }
   void set_mDx2(const Eigen::Matrix<double,3,3>& mDx2_) { mDx2 = mDx2_; }
   void set_mDx2(int i, int k, double value) { mDx2(i,k) = value; }
   void set_mDxbar2(const Eigen::Matrix<double,3,3>& mDxbar2_) { mDxbar2 = mDxbar2_; }
   void set_mDxbar2(int i, int k, double value) { mDxbar2(i,k) = value; }
   void set_mHp2(double mHp2_) { mHp2 = mHp2_; }
   void set_mHpbar2(double mHpbar2_) { mHpbar2 = mHpbar2_; }
   void set_mphi2(double mphi2_) { mphi2 = mphi2_; }
   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }
   void set_MassBp(double MassBp_) { MassBp = MassBp_; }

   const Eigen::Matrix<double,3,3>& get_TYd() const { return TYd; }
   double get_TYd(int i, int k) const { return TYd(i,k); }
   const Eigen::Matrix<double,3,2>& get_ThE() const { return ThE; }
   double get_ThE(int i, int k) const { return ThE(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYe() const { return TYe; }
   double get_TYe(int i, int k) const { return TYe(i,k); }
   double get_TSigmaL() const { return TSigmaL; }
   double get_TKappaPr() const { return TKappaPr; }
   double get_TSigmax() const { return TSigmax; }
   const Eigen::Matrix<double,3,3>& get_TgD() const { return TgD; }
   double get_TgD(int i, int k) const { return TgD(i,k); }
   const Eigen::Matrix<double,3,3>& get_TKappa() const { return TKappa; }
   double get_TKappa(int i, int k) const { return TKappa(i,k); }
   const Eigen::Matrix<double,2,2>& get_TLambda12() const { return TLambda12; }
   double get_TLambda12(int i, int k) const { return TLambda12(i,k); }
   double get_TLambdax() const { return TLambdax; }
   const Eigen::Matrix<double,3,2>& get_Tfu() const { return Tfu; }
   double get_Tfu(int i, int k) const { return Tfu(i,k); }
   const Eigen::Matrix<double,3,2>& get_Tfd() const { return Tfd; }
   double get_Tfd(int i, int k) const { return Tfd(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYu() const { return TYu; }
   double get_TYu(int i, int k) const { return TYu(i,k); }
   double get_BMuPr() const { return BMuPr; }
   double get_BMuPhi() const { return BMuPhi; }
   double get_LXiF() const { return LXiF; }
   const Eigen::Matrix<double,3,3>& get_mq2() const { return mq2; }
   double get_mq2(int i, int k) const { return mq2(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2() const { return ml2; }
   double get_ml2(int i, int k) const { return ml2(i,k); }
   double get_mHd2() const { return mHd2; }
   double get_mHu2() const { return mHu2; }
   const Eigen::Matrix<double,3,3>& get_md2() const { return md2; }
   double get_md2(int i, int k) const { return md2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2() const { return mu2; }
   double get_mu2(int i, int k) const { return mu2(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2() const { return me2; }
   double get_me2(int i, int k) const { return me2(i,k); }
   double get_ms2() const { return ms2; }
   double get_msbar2() const { return msbar2; }
   const Eigen::Matrix<double,2,2>& get_mH1I2() const { return mH1I2; }
   double get_mH1I2(int i, int k) const { return mH1I2(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH2I2() const { return mH2I2; }
   double get_mH2I2(int i, int k) const { return mH2I2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mSI2() const { return mSI2; }
   double get_mSI2(int i, int k) const { return mSI2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDx2() const { return mDx2; }
   double get_mDx2(int i, int k) const { return mDx2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2() const { return mDxbar2; }
   double get_mDxbar2(int i, int k) const { return mDxbar2(i,k); }
   double get_mHp2() const { return mHp2; }
   double get_mHpbar2() const { return mHpbar2; }
   double get_mphi2() const { return mphi2; }
   double get_MassB() const { return MassB; }
   double get_MassWB() const { return MassWB; }
   double get_MassG() const { return MassG; }
   double get_MassBp() const { return MassBp; }


protected:
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,2> ThE;
   Eigen::Matrix<double,3,3> TYe;
   double TSigmaL;
   double TKappaPr;
   double TSigmax;
   Eigen::Matrix<double,3,3> TgD;
   Eigen::Matrix<double,3,3> TKappa;
   Eigen::Matrix<double,2,2> TLambda12;
   double TLambdax;
   Eigen::Matrix<double,3,2> Tfu;
   Eigen::Matrix<double,3,2> Tfd;
   Eigen::Matrix<double,3,3> TYu;
   double BMuPr;
   double BMuPhi;
   double LXiF;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double ms2;
   double msbar2;
   Eigen::Matrix<double,2,2> mH1I2;
   Eigen::Matrix<double,2,2> mH2I2;
   Eigen::Matrix<double,3,3> mSI2;
   Eigen::Matrix<double,3,3> mDx2;
   Eigen::Matrix<double,3,3> mDxbar2;
   double mHp2;
   double mHpbar2;
   double mphi2;
   double MassB;
   double MassWB;
   double MassG;
   double MassBp;


private:
   static const int numberOfParameters = 249;

   struct Soft_traces {
      double traceAdjfdTfd;
      double traceAdjYdTYd;
      double traceAdjYeTYe;
      double tracefdAdjfd;
      double traceYdAdjYd;
      double traceYeAdjYe;
      double tracefuAdjfu;
      double traceYuAdjYu;
      double traceKappaAdjKappa;
      double traceLambda12AdjLambda12;
      double traceAdjfuTfu;
      double traceAdjYuTYu;
      double traceAdjKappaTKappa;
      double traceAdjLambda12TLambda12;
      double tracefdAdjfdTfdAdjfd;
      double tracefdAdjfdTfuAdjfu;
      double tracefuAdjfuTfdAdjfd;
      double tracehEAdjhETYeAdjYe;
      double traceYdAdjYdTYdAdjYd;
      double traceYdAdjYuTYuAdjYd;
      double traceYeAdjYeThEAdjhE;
      double traceYeAdjYeTYeAdjYe;
      double traceYuAdjYdTYdAdjYu;
      double traceAdjfdTfdconjLambda12TpLambda12;
      double traceAdjgDTpYdconjYdTgD;
      double traceAdjYdTYdconjgDTpgD;
      double traceAdjLambda12TpfdconjfdTLambda12;
      double traceAdjgDTgD;
      double traceAdjhEThE;
      double tracegDAdjgD;
      double tracehEAdjhE;
      double tracefdAdjfdfdAdjfd;
      double tracefdAdjfdfuAdjfu;
      double tracegDAdjgDTpYdconjYd;
      double tracehEAdjhEYeAdjYe;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceLambda12AdjLambda12Tpfdconjfd;
      double tracefuAdjhEThEAdjfu;
      double tracegDAdjgDTgDAdjgD;
      double tracegDAdjKappaTKappaAdjgD;
      double tracehEAdjfuTfuAdjhE;
      double tracehEAdjhEThEAdjhE;
      double tracehEAdjLambda12TLambda12AdjhE;
      double traceKappaAdjgDTgDAdjKappa;
      double traceLambda12AdjhEThEAdjLambda12;
      double traceAdjgDTpYuconjYuTgD;
      double traceAdjYuTYuconjgDTpgD;
      double tracefuAdjhEhEAdjfu;
      double tracegDAdjgDgDAdjgD;
      double tracegDAdjgDTpYuconjYu;
      double tracegDAdjKappaKappaAdjgD;
      double tracehEAdjhEhEAdjhE;
      double tracehEAdjLambda12Lambda12AdjhE;
      double tracefuAdjLambda12Lambda12Adjfu;
      double tracefuAdjLambda12TLambda12Adjfu;
      double traceKappaAdjKappaKappaAdjKappa;
      double traceKappaAdjKappaTKappaAdjKappa;
      double traceLambda12AdjfuTfuAdjLambda12;
      double traceLambda12AdjLambda12Lambda12AdjLambda12;
      double traceLambda12AdjLambda12TLambda12AdjLambda12;
      double tracefuAdjfufuAdjfu;
      double tracefuAdjfuTfuAdjfu;
      double traceYuAdjYuYuAdjYu;
      double traceYuAdjYuTYuAdjYu;
      double traceconjTgDTpgD;
      double traceconjTgDTpTgD;
      double traceconjThETphE;
      double traceconjThETpThE;
      double tracegDmDxbar2AdjgD;
      double tracegDAdjgDmq2;
      double tracehEAdjhEconjme2;
      double tracehEconjmH1I2AdjhE;
      double traceconjTfdTpTfd;
      double traceconjTYdTpTYd;
      double traceconjTYeTpTYe;
      double tracefdmH2I2Adjfd;
      double tracefdAdjfdconjmSI2;
      double tracemd2YdAdjYd;
      double traceme2YeAdjYe;
      double traceml2AdjYeYe;
      double tracemq2AdjYdYd;
      double traceconjTfdTpfd;
      double traceconjTYdTpYd;
      double traceconjTYeTpYe;
      double traceconjTfuTpTfu;
      double traceconjTYuTpTYu;
      double tracefumH1I2Adjfu;
      double tracefuAdjfuconjmSI2;
      double tracemq2AdjYuYu;
      double tracemu2YuAdjYu;
      double traceconjTfuTpfu;
      double traceconjTYuTpYu;
      double tracegDAdjgDconjmq2;
      double tracegDconjmDxbar2AdjgD;
      double tracehEmH1I2AdjhE;
      double tracehEAdjhEme2;
      double traceconjTKappaTpKappa;
      double traceconjTKappaTpTKappa;
      double traceconjTLambda12TpLambda12;
      double traceconjTLambda12TpTLambda12;
      double tracemH1I2AdjLambda12Lambda12;
      double traceKappaAdjKappaconjmDx2;
      double traceKappaconjmDxbar2AdjKappa;
      double traceLambda12AdjLambda12conjmH2I2;
      double tracefdAdjfdTfdAdjTfd;
      double tracefdAdjfdTfuAdjTfu;
      double tracefdAdjTfdTfdAdjfd;
      double tracefdAdjTfdTfuAdjfu;
      double tracefdconjTLambda12TpTLambda12Adjfd;
      double tracefuAdjfuTfdAdjTfd;
      double tracefuAdjTfuTfdAdjfd;
      double tracegDAdjgDTpTYdconjTYd;
      double tracehEAdjhETYeAdjTYe;
      double tracehEAdjThETYeAdjYe;
      double traceYdAdjYdTYdAdjTYd;
      double traceYdAdjYuTYuAdjTYd;
      double traceYdAdjTYdTYdAdjYd;
      double traceYdAdjTYuTYuAdjYd;
      double traceYdconjTgDTpTgDAdjYd;
      double traceYeAdjYeThEAdjThE;
      double traceYeAdjYeTYeAdjTYe;
      double traceYeAdjTYeThEAdjhE;
      double traceYeAdjTYeTYeAdjYe;
      double traceYuAdjYdTYdAdjTYu;
      double traceYuAdjTYdTYdAdjYu;
      double traceLambda12AdjLambda12TpTfdconjTfd;
      double traceAdjfdTfdconjTLambda12TpLambda12;
      double traceAdjgDTpYdconjTYdTgD;
      double traceAdjYdTYdconjTgDTpgD;
      double traceAdjLambda12TpfdconjTfdTLambda12;
      double tracefdmH2I2AdjfdfdAdjfd;
      double tracefdmH2I2AdjfdfuAdjfu;
      double tracefdAdjfdfdmH2I2Adjfd;
      double tracefdAdjfdfumH1I2Adjfu;
      double tracefdAdjfdfuAdjfuconjmSI2;
      double tracefdAdjfdconjmSI2fdAdjfd;
      double tracefdAdjfdconjmSI2fuAdjfu;
      double tracefdconjLambda12TpLambda12AdjfdconjmSI2;
      double tracegDAdjgDconjmq2TpYdconjYd;
      double tracegDAdjgDTpYdconjmd2conjYd;
      double tracegDAdjgDTpYdconjYdconjmq2;
      double tracegDconjmDxbar2AdjgDTpYdconjYd;
      double tracehEmH1I2AdjhEYeAdjYe;
      double tracehEAdjhEme2YeAdjYe;
      double tracehEAdjhEYeml2AdjYe;
      double tracehEAdjhEYeAdjYeme2;
      double tracemd2YdAdjYdYdAdjYd;
      double tracemd2YdAdjYuYuAdjYd;
      double traceme2YeAdjYeYeAdjYe;
      double tracemH1I2AdjLambda12TpfdconjfdLambda12;
      double traceml2AdjYeYeAdjYeYe;
      double tracemq2AdjYdYdAdjYdYd;
      double tracemq2AdjYdYdAdjYuYu;
      double tracemq2AdjYuYuAdjYdYd;
      double tracemu2YuAdjYdYdAdjYu;
      double traceLambda12AdjLambda12conjmH2I2Tpfdconjfd;
      double traceLambda12AdjLambda12TpfdconjfdconjmH2I2;
      double tracefuAdjfuTfuAdjTfu;
      double tracefuAdjhEThEAdjTfu;
      double tracefuAdjLambda12TLambda12AdjTfu;
      double tracefuAdjTfuTfuAdjfu;
      double tracefuAdjThEThEAdjfu;
      double tracefuAdjTLambda12TLambda12Adjfu;
      double tracegDAdjgDTpTYuconjTYu;
      double tracehEAdjfuTfuAdjThE;
      double tracehEAdjTfuTfuAdjhE;
      double traceYuAdjYuTYuAdjTYu;
      double traceYuAdjTYuTYuAdjYu;
      double traceYuconjTgDTpTgDAdjYu;
      double traceLambda12AdjfuTfuAdjTLambda12;
      double traceLambda12AdjTfuTfuAdjLambda12;
      double traceAdjgDTpYuconjTYuTgD;
      double traceAdjYuTYuconjTgDTpgD;
      double tracefumH1I2AdjfufuAdjfu;
      double tracefumH1I2AdjhEhEAdjfu;
      double tracefumH1I2AdjLambda12Lambda12Adjfu;
      double tracefuAdjfufumH1I2Adjfu;
      double tracefuAdjfuconjmSI2fuAdjfu;
      double tracefuAdjhEhEmH1I2Adjfu;
      double tracefuAdjhEhEAdjfuconjmSI2;
      double tracefuAdjhEme2hEAdjfu;
      double tracefuAdjLambda12Lambda12mH1I2Adjfu;
      double tracefuAdjLambda12Lambda12AdjfuconjmSI2;
      double tracefuAdjLambda12conjmH2I2Lambda12Adjfu;
      double tracegDAdjgDconjmq2TpYuconjYu;
      double tracegDAdjgDTpYuconjmu2conjYu;
      double tracegDAdjgDTpYuconjYuconjmq2;
      double tracegDconjmDxbar2AdjgDTpYuconjYu;
      double tracemq2AdjYuYuAdjYuYu;
      double tracemu2YuAdjYuYuAdjYu;
      double tracegDAdjKappaTKappaAdjTgD;
      double tracegDAdjTKappaTKappaAdjgD;
      double tracehEAdjLambda12TLambda12AdjThE;
      double tracehEAdjTLambda12TLambda12AdjhE;
      double traceKappaAdjgDTgDAdjTKappa;
      double traceKappaAdjKappaTKappaAdjTKappa;
      double traceKappaAdjTgDTgDAdjKappa;
      double traceKappaAdjTKappaTKappaAdjKappa;
      double traceLambda12AdjhEThEAdjTLambda12;
      double traceLambda12AdjLambda12TLambda12AdjTLambda12;
      double traceLambda12AdjThEThEAdjLambda12;
      double traceLambda12AdjTLambda12TLambda12AdjLambda12;
      double tracegDAdjKappaKappaAdjgDconjmq2;
      double tracegDAdjKappaKappaconjmDxbar2AdjgD;
      double tracegDAdjKappaconjmDx2KappaAdjgD;
      double tracegDconjmDxbar2AdjKappaKappaAdjgD;
      double tracehEmH1I2AdjLambda12Lambda12AdjhE;
      double tracehEAdjLambda12Lambda12mH1I2AdjhE;
      double tracehEAdjLambda12Lambda12AdjhEme2;
      double tracehEAdjLambda12conjmH2I2Lambda12AdjhE;
      double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
      double traceKappaAdjKappaKappaAdjKappaconjmDx2;
      double traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
      double traceKappaAdjKappaconjmDx2KappaAdjKappa;
      double traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
      double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
      double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
      double tracegDAdjgDTgDAdjTgD;
      double tracegDAdjTgDTgDAdjgD;
      double tracehEAdjhEThEAdjThE;
      double tracehEAdjThEThEAdjhE;
      double tracegDAdjgDgDconjmDxbar2AdjgD;
      double tracegDAdjgDconjmq2gDAdjgD;
      double tracegDconjmDxbar2AdjgDgDAdjgD;
      double tracehEmH1I2AdjhEhEAdjhE;
      double tracehEAdjhEhEmH1I2AdjhE;
      double tracehEAdjhEhEAdjhEme2;
      double tracehEAdjhEme2hEAdjhE;
      double Tr11;
      double Tr14;
      double Tr2U111;
      double Tr2U114;
      double Tr31;
      double Tr22;
      double Tr23;
      double Tr2U141;
      double Tr2U144;
      double Tr34;

   };
   void calc_soft_traces(Soft_traces&) const;

   Eigen::Matrix<double,3,3> calc_beta_TYd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_ThE_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_ThE_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmaL_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmaL_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappaPr_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappaPr_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TgD_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TgD_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TKappa_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TKappa_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_TLambda12_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_TLambda12_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_Tfu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_Tfu_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_Tfd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,2> calc_beta_Tfd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPr_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPr_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPhi_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPhi_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_LXiF_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_LXiF_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_msbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_msbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_mH1I2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_mH1I2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_mH2I2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_mH2I2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mSI2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mSI2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDx2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDx2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDxbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDxbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHp2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHp2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHpbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHpbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mphi2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mphi2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassBp_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassBp_two_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const CSE6SSM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
