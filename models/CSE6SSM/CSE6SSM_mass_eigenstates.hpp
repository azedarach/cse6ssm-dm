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

// File generated at Wed 3 Jun 2015 23:53:01

/**
 * @file CSE6SSM_mass_eigenstates.hpp
 * @brief contains class for model with routines needed for
 *        determining the pole masses and mixings
 *
 * This file was generated at Wed 3 Jun 2015 23:53:01 with FlexibleSUSY
 * 1.1.0 (git commit: v1.1.0) and SARAH 4.5.6 .
 */

#ifndef CSE6SSM_MASS_EIGENSTATES_H
#define CSE6SSM_MASS_EIGENSTATES_H

#include "CSE6SSM_soft_parameters.hpp"
#include "CSE6SSM_physical.hpp"
#include "CSE6SSM_info.hpp"
#include "two_loop_corrections.hpp"
#include "problems.hpp"
#include "config.h"

#ifdef ENABLE_THREADS
#include <mutex>
#endif

#include <iosfwd>

#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class CSE6SSM_mass_eigenstates
 * @brief model class with routines for determing masses and mixings
 */
class CSE6SSM_mass_eigenstates : public CSE6SSM_soft_parameters {
public:
   explicit CSE6SSM_mass_eigenstates();
   virtual ~CSE6SSM_mass_eigenstates();

   /// number of EWSB equations
   static const std::size_t number_of_tadpole_equations = 5;

   void calculate_DRbar_masses();
   void calculate_DRbar_parameters();
   void calculate_pole_masses();
   virtual void clear();
   void clear_DRbar_parameters();
   void do_calculate_sm_pole_masses(bool);
   bool do_calculate_sm_pole_masses() const;
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_two_loop_corrections(const Two_loop_corrections&);
   void set_number_of_mass_iterations(std::size_t);
   void set_diagonalization_precision(double);
   double get_diagonalization_precision() const;
   void set_pole_mass_loop_order(unsigned);
   void set_physical(const CSE6SSM_physical&);
   const CSE6SSM_physical& get_physical() const;
   CSE6SSM_physical& get_physical();
   const CSE6SSM_physical& get_drbar_masses() const;
   CSE6SSM_physical& get_drbar_masses();
   const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& get_problems() const;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& get_problems();
   const Eigen::Array<double,5,1> get_ewsb_tree_level_soft_masses();

   virtual void print(std::ostream&) const;

   double get_MVG() const { return MVG; }
   double get_MGlu() const { return MGlu; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_MChaP() const { return MChaP; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }
   double get_MVZp() const { return MVZp; }
   const Eigen::Array<double,6,1>& get_MSd() const { return MSd; }
   double get_MSd(int i) const { return MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv() const { return MSv; }
   double get_MSv(int i) const { return MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu() const { return MSu; }
   double get_MSu(int i) const { return MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe() const { return MSe; }
   double get_MSe(int i) const { return MSe(i); }
   const Eigen::Array<double,6,1>& get_MSDX() const { return MSDX; }
   double get_MSDX(int i) const { return MSDX(i); }
   const Eigen::Array<double,5,1>& get_Mhh() const { return Mhh; }
   double get_Mhh(int i) const { return Mhh(i); }
   const Eigen::Array<double,5,1>& get_MAh() const { return MAh; }
   double get_MAh(int i) const { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm() const { return MHpm; }
   double get_MHpm(int i) const { return MHpm(i); }
   const Eigen::Array<double,8,1>& get_MChi() const { return MChi; }
   double get_MChi(int i) const { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha() const { return MCha; }
   double get_MCha(int i) const { return MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFDX() const { return MFDX; }
   double get_MFDX(int i) const { return MFDX(i); }
   const Eigen::Array<double,7,1>& get_MSHI0() const { return MSHI0; }
   double get_MSHI0(int i) const { return MSHI0(i); }
   const Eigen::Array<double,4,1>& get_MSHIPM() const { return MSHIPM; }
   double get_MSHIPM(int i) const { return MSHIPM(i); }
   const Eigen::Array<double,2,1>& get_MChaI() const { return MChaI; }
   double get_MChaI(int i) const { return MChaI(i); }
   const Eigen::Array<double,7,1>& get_MChiI() const { return MChiI; }
   double get_MChiI(int i) const { return MChiI(i); }
   const Eigen::Array<double,2,1>& get_MSHp0() const { return MSHp0; }
   double get_MSHp0(int i) const { return MSHp0(i); }
   const Eigen::Array<double,2,1>& get_MSHpp() const { return MSHpp; }
   double get_MSHpp(int i) const { return MSHpp(i); }
   const Eigen::Array<double,2,1>& get_MChiP() const { return MChiP; }
   double get_MChiP(int i) const { return MChiP(i); }
   double get_MVWm() const { return MVWm; }

   
   Eigen::Array<double,1,1> get_MChargedHiggs() const;

   Eigen::Array<double,3,1> get_MPseudoscalarHiggs() const;

   const Eigen::Matrix<double,6,6>& get_ZD() const { return ZD; }
   double get_ZD(int i, int k) const { return ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV() const { return ZV; }
   double get_ZV(int i, int k) const { return ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU() const { return ZU; }
   double get_ZU(int i, int k) const { return ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE() const { return ZE; }
   double get_ZE(int i, int k) const { return ZE(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZDX() const { return ZDX; }
   double get_ZDX(int i, int k) const { return ZDX(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZH() const { return ZH; }
   double get_ZH(int i, int k) const { return ZH(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZA() const { return ZA; }
   double get_ZA(int i, int k) const { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const { return ZP; }
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,8,8>& get_ZN() const { return ZN; }
   const std::complex<double>& get_ZN(int i, int k) const { return ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM() const { return UM; }
   const std::complex<double>& get_UM(int i, int k) const { return UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP() const { return UP; }
   const std::complex<double>& get_UP(int i, int k) const { return UP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL() const { return ZEL; }
   const std::complex<double>& get_ZEL(int i, int k) const { return ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER() const { return ZER; }
   const std::complex<double>& get_ZER(int i, int k) const { return ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL() const { return ZDL; }
   const std::complex<double>& get_ZDL(int i, int k) const { return ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR() const { return ZDR; }
   const std::complex<double>& get_ZDR(int i, int k) const { return ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL() const { return ZUL; }
   const std::complex<double>& get_ZUL(int i, int k) const { return ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR() const { return ZUR; }
   const std::complex<double>& get_ZUR(int i, int k) const { return ZUR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXL() const { return ZDXL; }
   const std::complex<double>& get_ZDXL(int i, int k) const { return ZDXL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXR() const { return ZDXR; }
   const std::complex<double>& get_ZDXR(int i, int k) const { return ZDXR(i,k); }
   const Eigen::Matrix<double,7,7>& get_UHI0() const { return UHI0; }
   double get_UHI0(int i, int k) const { return UHI0(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHIPM() const { return UHIPM; }
   double get_UHIPM(int i, int k) const { return UHIPM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZMI() const { return ZMI; }
   const std::complex<double>& get_ZMI(int i, int k) const { return ZMI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZPI() const { return ZPI; }
   const std::complex<double>& get_ZPI(int i, int k) const { return ZPI(i,k); }
   const Eigen::Matrix<std::complex<double>,7,7>& get_ZNI() const { return ZNI; }
   const std::complex<double>& get_ZNI(int i, int k) const { return ZNI(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHp0() const { return UHp0; }
   double get_UHp0(int i, int k) const { return UHp0(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHpp() const { return UHpp; }
   double get_UHpp(int i, int k) const { return UHpp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZNp() const { return ZNp; }
   const std::complex<double>& get_ZNp(int i, int k) const { return ZNp(i,k); }

   void set_PhaseGlu(std::complex<double> PhaseGlu_) { PhaseGlu = PhaseGlu_; }
   std::complex<double> get_PhaseGlu() const { return PhaseGlu; }

   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_Glu() const;
   void calculate_MGlu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   double get_mass_matrix_ChaP() const;
   void calculate_MChaP();
   double get_mass_matrix_VP() const;
   void calculate_MVP();
   double get_mass_matrix_VZ() const;
   void calculate_MVZ();
   double get_mass_matrix_VZp() const;
   void calculate_MVZp();
   Eigen::Matrix<double,6,6> get_mass_matrix_Sd() const;
   void calculate_MSd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Sv() const;
   void calculate_MSv();
   Eigen::Matrix<double,6,6> get_mass_matrix_Su() const;
   void calculate_MSu();
   Eigen::Matrix<double,6,6> get_mass_matrix_Se() const;
   void calculate_MSe();
   Eigen::Matrix<double,6,6> get_mass_matrix_SDX() const;
   void calculate_MSDX();
   Eigen::Matrix<double,5,5> get_mass_matrix_hh() const;
   Eigen::Matrix<double,5,5> get_rotated_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,5,5> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Hpm() const;
   void calculate_MHpm();
   Eigen::Matrix<double,8,8> get_mass_matrix_Chi() const;
   void calculate_MChi();
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const;
   void calculate_MCha();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<double,3,3> get_mass_matrix_FDX() const;
   void calculate_MFDX();
   Eigen::Matrix<double,7,7> get_mass_matrix_SHI0() const;
   void calculate_MSHI0();
   Eigen::Matrix<double,4,4> get_mass_matrix_SHIPM() const;
   void calculate_MSHIPM();
   Eigen::Matrix<double,2,2> get_mass_matrix_ChaI() const;
   void calculate_MChaI();
   Eigen::Matrix<double,7,7> get_mass_matrix_ChiI() const;
   void calculate_MChiI();
   Eigen::Matrix<double,2,2> get_mass_matrix_SHp0() const;
   void calculate_MSHp0();
   Eigen::Matrix<double,2,2> get_mass_matrix_SHpp() const;
   void calculate_MSHpp();
   Eigen::Matrix<double,2,2> get_mass_matrix_ChiP() const;
   void calculate_MChiP();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;
   double get_ewsb_eq_hh_3() const;
   double get_ewsb_eq_hh_4() const;
   double get_ewsb_eq_hh_5() const;

   double get_tadpole_eq_hh_1(unsigned ewsb_loop_order) const;
   double get_tadpole_eq_hh_2(unsigned ewsb_loop_order) const;
   double get_tadpole_eq_hh_3(unsigned ewsb_loop_order) const;
   double get_tadpole_eq_hh_4(unsigned ewsb_loop_order) const;
   double get_tadpole_eq_hh_5(unsigned ewsb_loop_order) const;

   std::complex<double> CpUSdconjUSdVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSdconjUSdVZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSdconjUSdconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdconjSHp0SDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFDXChiPPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSdFDXChiPPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSdGluFdPL(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSdVGSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVPSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVZSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVZpSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVWmSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSvconjUSvVZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvbarChaFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSvFvChiPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvVZSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSvVZpSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSvconjVWmSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSuconjUSuVZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSuconjUSuconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarChaFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuconjSHppSDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarChaPFDXPR(unsigned gO2, unsigned gI2) const;
   double CpconjUSubarChaPFDXPL(unsigned , unsigned ) const;
   std::complex<double> CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSuGluFuPL(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSuconjVWmSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVGSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVPSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVZSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVZpSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSeconjUSeVZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSeconjUSeconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSeChiPChaIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSeChiPChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSHp0SHIPM(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSHppSHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSeFvChaPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSeFvChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeconjSHI0SHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSeChiIChaPPR(unsigned , unsigned ) const;
   std::complex<double> CpconjUSeChiIChaPPL(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpconjUSeVWmSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVPSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVZSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVZpSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSDXconjUSDXVZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSDXconjUSDXconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXSHp0Sd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSDXFdChiPPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSDXFdChiPPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXFDXChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXFDXChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpconjUSDXFuChaPPR(unsigned , unsigned ) const;
   std::complex<double> CpconjUSDXFuChaPPL(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpUSDXconjUSDXconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXSuSHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXSDXAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXSDXhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSDXGluFDXPR(unsigned gO2, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSDXGluFDXPL(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSDXVGSDX(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSDXVPSDX(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSDXVZSDX(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSDXVZpSDX(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhVZVZ(unsigned gO2) const;
   std::complex<double> CpUhhVZpVZ(unsigned gO2) const;
   std::complex<double> CpUhhVZpVZp(unsigned gO2) const;
   std::complex<double> CpUhhconjVWmVWm(unsigned gO2) const;
   std::complex<double> CpUhhbargWmgWm(unsigned gO1) const;
   std::complex<double> CpUhhbargWmCgWmC(unsigned gO1) const;
   std::complex<double> CpUhhbargZgZ(unsigned gO1) const;
   std::complex<double> CpUhhbargZpgZ(unsigned gO1) const;
   std::complex<double> CpUhhbargZpgZp(unsigned gO1) const;
   std::complex<double> CpUhhUhhVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhVZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSHp0SHp0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSHppSHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChaIChaIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChaIChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiPChiPPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiPChiPPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFDXFDXPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFDXFDXPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSHIPMSHIPM(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSDXSDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSHI0SHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhSHI0SHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiIChiIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiIChiIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjVWmHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhVZAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhVZpAh(unsigned gO2, unsigned gI2) const;
   double CpUhhbarChaPChaPPR(unsigned gO2) const;
   double CpUhhbarChaPChaPPL(unsigned gO1) const;
   std::complex<double> CpUAhbargWmgWm(unsigned gO1) const;
   std::complex<double> CpUAhbargWmCgWmC(unsigned gO1) const;
   std::complex<double> CpUAhUAhVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhVZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSHp0SHp0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSHppSHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChaIChaIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChaIChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiPChiPPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiPChiPPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFDXFDXPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFDXFDXPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSHIPMSHIPM(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSDXSDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSHI0SHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhSHI0SHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiIChiIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiIChiIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjVWmHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhVZhh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhVZphh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhbarChaPChaPPR(unsigned gO2) const;
   std::complex<double> CpUAhbarChaPChaPPL(unsigned gO1) const;
   std::complex<double> CpconjUHpmVWmVP(unsigned gO2) const;
   std::complex<double> CpconjUHpmVZVWm(unsigned gO2) const;
   std::complex<double> CpconjUHpmVZpVWm(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargWmCgZ(unsigned gO1) const;
   std::complex<double> CpUHpmgWmCbargZ(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargWmCgZp(unsigned gO1) const;
   std::complex<double> CpUHpmgWmCbargZp(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargZgWm(unsigned gO1) const;
   std::complex<double> CpUHpmgZbargWm(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargZpgWm(unsigned gO1) const;
   std::complex<double> CpUHpmgZpbargWm(unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmVZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSHp0SHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUHpmbarFvFePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmSHIPMSHI0(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmChiIChaIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmChiIChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSHI0SHIPM(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmVPHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVZHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVZpHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0VZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHI0conjUSHI0VZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSHI0conjUSHI0conjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHI0conjUSHI0conjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0barChaChaIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0barChaChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0barChaIChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHI0barChaIChaPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSHI0conjHpmSHIPM(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0conjSHppSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0conjSHIPMHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0AhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0hhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0conjSeSHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHI0conjUSHI0conjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0conjSHI0Ah(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0conjSHI0hh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0SHI0Ah(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0SHI0hh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0ChiIChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0ChiIChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHI0barChaPFePR(unsigned gO2, unsigned gI2) const;
   double CpconjUSHI0barChaPFePL(unsigned , unsigned ) const;
   std::complex<double> CpconjUSHI0conjVWmSHIPM(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHI0VZSHI0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHI0VZpSHI0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHIPMconjUSHIPMVZpVZp(unsigned gO1, unsigned gO2) const;
   double CpUSHIPMconjUSHIPMconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMconjSHp0Se(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMFeChiPPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHIPMFeChiPPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMSHIPMAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMSHIPMhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHIPMconjUSHIPMconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMconjSHI0Hpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMSHI0Hpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMChiIChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMChiIChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMChiChaIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMChiChaIPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMVPSHIPM(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMVZSHIPM(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMVZpSHIPM(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHIPMVWmSHI0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0VZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHp0conjUSHp0VZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHp0conjUSHp0conjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHp0conjUSHp0conjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0conjHpmSHpp(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0barChaIFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHp0barChaIFePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSHp0SHp0Ah(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0SHp0hh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0ChiPChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0ChiPChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0barChaChaPPR(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpconjUSHp0barChaChaPPL(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpUSHp0conjUSHp0conjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0barFdFDXPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHp0barFdFDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUSHp0conjUSHp0conjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0conjSHIPMSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0AhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0hhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0conjSdSDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHp0conjUSHp0conjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHp0VZSHp0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHp0VZpSHp0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHp0conjVWmSHpp(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHppconjUSHppVZpVZp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHppconjUSHppconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSHppconjUSHppconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSHp0SHp0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSHppSHpp(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppSHp0Hpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppChiPChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppChiPChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppSHppAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppSHpphh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppbarFuFDXPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHppbarFuFDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSHppFeChiIPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSHppFeChiIPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUSHppconjUSHppconjSHIPMSHIPM(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHpphhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSDXSDX(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppconjSuSDX(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSHppconjUSHppconjSHI0SHI0(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppconjSHI0Se(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppSHI0Se(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSHppChiChaPPR(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpconjUSHppChiChaPPL(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpconjUSHppVWmSHp0(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHppVPSHpp(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHppVZSHpp(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSHppVZpSHpp(unsigned gO2, unsigned gI2) const;
   double CpVZbargWmgWm() const;
   double CpVZbargWmCgWmC() const;
   double CpVZconjVWmVWm() const;
   double CpVZbarChaPChaPPL() const;
   double CpVZbarChaPChaPPR() const;
   std::complex<double> CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSHp0SHp0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSHppSHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const;
   double CpVZconjSHp0SHp0(unsigned gI1, unsigned gI2) const;
   double CpVZconjSHppSHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarChaChaPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarChaChaPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarChaIChaIPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarChaIChaIPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiPChiPPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiPChiPPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFDXFDXPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFDXFDXPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFeFePL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFeFePR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFvFvPR(unsigned , unsigned ) const;
   std::complex<double> CpVZVZconjSHIPMSHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSHIPMSHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZhhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZhhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSDXSDX(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSDXSDX(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSHI0SHI0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSHI0SHI0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiIChiIPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiIChiIPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiChiPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZChiChiPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjVWmHpm(unsigned gI2) const;
   std::complex<double> CpVZVZhh(unsigned gI2) const;
   std::complex<double> CpVZVZphh(unsigned gI2) const;
   double CpVZVZconjVWmVWm1() const;
   double CpVZVZconjVWmVWm2() const;
   double CpVZVZconjVWmVWm3() const;
   double CpVZpbargWmgWm() const;
   double CpVZpbargWmCgWmC() const;
   double CpVZpconjVWmVWm() const;
   double CpVZpbarChaPChaPPL() const;
   double CpVZpbarChaPChaPPR() const;
   std::complex<double> CpVZpVZpconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSHp0SHp0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSHppSHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjHpmHpm(unsigned gI1, unsigned gI2) const;
   double CpVZpconjSHp0SHp0(unsigned gI1, unsigned gI2) const;
   double CpVZpconjSHppSHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpbarChaChaPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpbarChaChaPR(unsigned gI1, unsigned gI2) const;
   double CpVZpbarChaIChaIPL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarChaIChaIPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiPChiPPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiPChiPPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZpconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFdFdPL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFdFdPR(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFDXFDXPL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFDXFDXPR(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFeFePL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFeFePR(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFuFuPL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFuFuPR(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFvFvPL(unsigned gI1, unsigned gI2) const;
   double CpVZpbarFvFvPR(unsigned , unsigned ) const;
   std::complex<double> CpVZpVZpconjSHIPMSHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSHIPMSHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZphhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZphhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSDXSDX(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSDXSDX(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpVZpconjSHI0SHI0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjSHI0SHI0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiIChiIPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiIChiIPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiChiPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpChiChiPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZpconjVWmHpm(unsigned gI2) const;
   std::complex<double> CpVZpVZhh(unsigned gI2) const;
   std::complex<double> CpVZpVZphh(unsigned gI2) const;
   double CpVZpVZpconjVWmVWm1() const;
   double CpVZpVZpconjVWmVWm2() const;
   double CpVZpVZpconjVWmVWm3() const;
   double CpconjVWmbargPgWm() const;
   double CpconjVWmbargWmCgP() const;
   double CpconjVWmbargWmCgZ() const;
   double CpconjVWmbargWmCgZp() const;
   double CpconjVWmbargZgWm() const;
   double CpconjVWmbargZpgWm() const;
   double CpconjVWmVWmVP() const;
   double CpconjVWmVZVWm() const;
   double CpconjVWmVZpVWm() const;
   std::complex<double> CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSHp0SHp0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSHppSHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmconjSHp0SHpp(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmChiPChaPPL(unsigned gI1) const;
   std::complex<double> CpconjVWmChiPChaPPR(unsigned gI1) const;
   double CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const;
   double CpconjVWmbarFuFdPR(unsigned , unsigned ) const;
   std::complex<double> CpconjVWmbarFvFePL(unsigned gI1, unsigned gI2) const;
   double CpconjVWmbarFvFePR(unsigned , unsigned ) const;
   std::complex<double> CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const;
   double CpVWmconjVWmconjSHIPMSHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSHI0SHI0(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmChiIChaIPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmChiIChaIPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmconjSHI0SHIPM(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmChiChaPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmChiChaPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmVPHpm(unsigned gI2) const;
   std::complex<double> CpconjVWmVZHpm(unsigned gI2) const;
   std::complex<double> CpconjVWmVZpHpm(unsigned gI2) const;
   std::complex<double> CpconjVWmVWmhh(unsigned gI2) const;
   double CpVWmconjVWmVPVP1() const;
   double CpVWmconjVWmVPVP2() const;
   double CpVWmconjVWmVPVP3() const;
   double CpVWmconjVWmVZVZ1() const;
   double CpVWmconjVWmVZVZ2() const;
   double CpVWmconjVWmVZVZ3() const;
   double CpVWmconjVWmVZpVZp1() const;
   double CpVWmconjVWmVZpVZp2() const;
   double CpVWmconjVWmVZpVZp3() const;
   double CpVWmconjVWmconjVWmVWm1() const;
   double CpVWmconjVWmconjVWmVWm2() const;
   double CpVWmconjVWmconjVWmVWm3() const;
   std::complex<double> CpUChiconjHpmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjHpmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHp0ChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHp0ChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHppChaPPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpUChiconjSHppChaPPR(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpUChiconjSvFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUChiconjSHIPMChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHIPMChaIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSDXFDXPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSDXFDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHI0ChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjSHI0ChiIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiconjVWmChaPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiconjVWmChaPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiVZChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiVZChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiVZpChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiVZpChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaSHppChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaSHppChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSHp0ChaPPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpbarUChaconjSHp0ChaPPR(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUChabarFvSePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChabarFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaSHIPMChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaSHIPMChiIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSHI0ChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaconjSHI0ChaIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaSHI0ChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUChaSHI0ChaIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChaVPChaPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaVPChaPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaVZChaPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaVZChaPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaVZpChaPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaVZpChaPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaVWmChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaVWmChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFeSHp0ChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFeSHp0ChaIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFeHpmFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFeSHppChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFeSHppChiIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSHIPMChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFeSHIPMChiPPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSHI0ChaPPL(unsigned gO2, unsigned gI1) const;
   double CpbarUFeSHI0ChaPPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const;
   double CpbarUFeVWmFvPR(unsigned , unsigned ) const;
   double CpbarUFeVWmFvPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFeVZpFePR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFeVZpFePL(unsigned gO1, unsigned gI2) const;
   double CpbarUFdconjSHp0FDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFdconjSHp0FDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUFdSDXChiPPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFdSDXChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFdVWmFuPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFdVWmFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdVZpFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVZpFdPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUFuconjSHppFDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFuconjSHppFDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuVZpFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVZpFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFuconjVWmFdPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFuconjVWmFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFubarChaPSDXPL(unsigned , unsigned ) const;
   std::complex<double> CpbarUFubarChaPSDXPR(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFDXSHp0FdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFDXSHp0FdPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFDXSHppFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFDXSHppFuPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFDXFDXAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXFDXAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXhhFDXPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXhhFDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXSdChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFDXSdChiPPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFDXSDXChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXSDXChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFDXSDXGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFDXSDXGluPR(unsigned gO1, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarUFDXSuChaPPL(unsigned gO2, unsigned gI1) const;
   double CpbarUFDXSuChaPPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFDXVGFDXPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFDXVGFDXPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFDXVPFDXPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFDXVPFDXPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFDXVZFDXPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFDXVZFDXPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFDXVZpFDXPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFDXVZpFDXPL(unsigned gO1, unsigned gI2) const;
   double CpbarUChaIconjSHp0FePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChaIconjSHp0FePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIChaIAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIChaIAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIHpmChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIHpmChiIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaISHIPMChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaISHIPMChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIhhChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIhhChaIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUChaISeChiPPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChaISeChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUChaIconjSHI0ChaPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChaIconjSHI0ChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaISHI0ChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaISHI0ChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChaIVPChaIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaIVPChaIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaIVZChaIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaIVZChaIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaIVZpChaIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaIVZpChaIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChaIVWmChiIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChaIVWmChiIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiIconjHpmChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjHpmChaIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpUChiIconjSHppFePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUChiIconjSHppFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjSHIPMChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjSHIPMChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIhhChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIhhChiIPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjSeChaPPL(unsigned gO2, unsigned gI1) const;
   double CpUChiIconjSeChaPPR(unsigned , unsigned ) const;
   std::complex<double> CpUChiIChiIAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIChiIAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjSHI0ChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjSHI0ChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiIconjVWmChaIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiIconjVWmChaIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiIVZChiIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiIVZChiIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiIVZpChiIPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiIVZpChiIPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSHppChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSHppChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPChiPAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPChiPAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSHp0ChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSHp0ChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpUChiPconjSHIPMFePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUChiPconjSHIPMFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPhhChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPhhChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSeChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpUChiPconjSeChaIPR(unsigned , unsigned , unsigned ) const;
   double CpUChiPconjSdFDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUChiPconjSdFDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUChiPconjSDXFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpUChiPconjSDXFdPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpUChiPVZChiPPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiPVZChiPPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpUChiPVZpChiPPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUChiPVZpChiPPL(unsigned gO1, unsigned gI2) const;
   double CpUChiPconjVWmChaPPR(unsigned gO2) const;
   double CpUChiPconjVWmChaPPL(unsigned gO1) const;
   std::complex<double> CpGluconjSdFdPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluconjSDXFDXPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluconjSDXFDXPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpGluVGGluPR() const;
   std::complex<double> CpGluVGGluPL() const;
   std::complex<double> CpbarChaPSHp0ChaPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarChaPSHp0ChaPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarChaPSHppChiPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarChaPSHppChiPR(unsigned gI1, unsigned gI2) const;
   double CpbarChaPbarFuSDXPL(unsigned , unsigned ) const;
   std::complex<double> CpbarChaPbarFuSDXPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarChaPhhChaPPL(unsigned gI1) const;
   std::complex<double> CpbarChaPhhChaPPR(unsigned gI1) const;
   double CpbarChaPconjSuFDXPL(unsigned , unsigned ) const;
   std::complex<double> CpbarChaPconjSuFDXPR(unsigned gI1, unsigned gI2) const;
   double CpbarChaPSeChiIPL(unsigned , unsigned ) const;
   std::complex<double> CpbarChaPSeChiIPR(unsigned gI1, unsigned gI2) const;
   double CpbarChaPconjSHI0FePL(unsigned , unsigned ) const;
   std::complex<double> CpbarChaPconjSHI0FePR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarChaPVWmChiPPR(unsigned gI2) const;
   std::complex<double> CpbarChaPVWmChiPPL(unsigned gI2) const;
   std::complex<double> CpbarChaPChaPAhPL(unsigned gI2) const;
   std::complex<double> CpbarChaPChaPAhPR(unsigned gI2) const;
   double CpbarChaPVPChaPPR() const;
   double CpbarChaPVPChaPPL() const;
   double CpbarChaPVZChaPPR() const;
   double CpbarChaPVZChaPPL() const;
   double CpbarChaPVZpChaPPR() const;
   double CpbarChaPVZpChaPPL() const;
   double CpconjVWmbarVWmVZp() const;
   double CpconjVWmbarVZpVWm() const;
   std::complex<double> CpbarFeSHp0ChaIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarFeSHp0ChaIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFeSHppChiIPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarFeSHppChiIPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSHIPMChiPPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarFeSHIPMChiPPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSHI0ChaPPL(unsigned gO2, unsigned gI1) const;
   double CpbarFeSHI0ChaPPR(unsigned , unsigned ) const;
   double CpbarFeVWmFvPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFeVWmFvPL(unsigned gO1, unsigned gI2) const;
   double CpbarFeVZFePR(unsigned gO2, unsigned gI2) const;
   double CpbarFeVZFePL(unsigned gO1, unsigned gI2) const;
   double CpbarFeVZpFePR(unsigned gO2, unsigned gI2) const;
   double CpbarFeVZpFePL(unsigned gO1, unsigned gI2) const;
   double CpbarFdconjSHp0FDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFdconjSHp0FDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarFdSDXChiPPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFdSDXChiPPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const;
   double CpbarFdVWmFuPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFdVWmFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const;
   double CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarFdVZpFdPR(unsigned gO2, unsigned gI2) const;
   double CpbarFdVZpFdPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarFuconjSHppFDXPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFuconjSHppFDXPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   std::complex<double> CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const;
   double CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const;
   double CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const;
   double CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFuVZpFuPR(unsigned gO2, unsigned gI2) const;
   double CpbarFuVZpFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFuconjVWmFdPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFuconjVWmFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarFubarChaPSDXPL(unsigned , unsigned ) const;
   std::complex<double> CpbarFubarChaPSDXPR(unsigned gO1, unsigned gI2) const;
   std::complex<double> self_energy_Sd(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Sv(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Su(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Se(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SDX(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_hh(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Ah(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SHI0(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SHIPM(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SHp0(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SHpp(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_VZ(double p ) const;
   std::complex<double> self_energy_VZp(double p ) const;
   std::complex<double> self_energy_VWm(double p ) const;
   std::complex<double> self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_FDX_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_FDX_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_FDX_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChaI_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChaI_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChaI_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiI_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiI_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiI_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiP_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiP_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_ChiP_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Glu_1(double p ) const;
   std::complex<double> self_energy_Glu_PR(double p ) const;
   std::complex<double> self_energy_Glu_PL(double p ) const;
   std::complex<double> self_energy_ChaP_1(double p ) const;
   std::complex<double> self_energy_ChaP_PR(double p ) const;
   std::complex<double> self_energy_ChaP_PL(double p ) const;
   std::complex<double> self_energy_VZ_heavy(double p ) const;
   std::complex<double> self_energy_VWm_heavy(double p ) const;
   std::complex<double> self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> tadpole_hh(unsigned gO1) const;

   void calculate_MSu_3rd_generation(double&, double&, double&) const;
   void calculate_MSd_3rd_generation(double&, double&, double&) const;
   void calculate_MSv_3rd_generation(double&, double&, double&) const;
   void calculate_MSe_3rd_generation(double&, double&, double&) const;
 
   void self_energy_hh_2loop(double result[6]) const;
   void self_energy_Ah_2loop(double result[6]) const;
 
   void tadpole_hh_2loop(double result[3]) const;

   void calculate_MVG_pole();
   void calculate_MGlu_pole();
   void calculate_MFv_pole();
   void calculate_MChaP_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MVZp_pole();
   void calculate_MSd_pole();
   void calculate_MSv_pole();
   void calculate_MSu_pole();
   void calculate_MSe_pole();
   void calculate_MSDX_pole();
   void calculate_Mhh_pole();
   void calculate_MAh_pole();
   void calculate_MHpm_pole();
   void calculate_MChi_pole();
   void calculate_MCha_pole();
   void calculate_MFe_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFDX_pole();
   void calculate_MSHI0_pole();
   void calculate_MSHIPM_pole();
   void calculate_MChaI_pole();
   void calculate_MChiI_pole();
   void calculate_MSHp0_pole();
   void calculate_MSHpp_pole();
   void calculate_MChiP_pole();
   void calculate_MVWm_pole();
   double calculate_MVWm_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFv_DRbar(double, int) const;
   double calculate_MVP_DRbar(double);
   double calculate_MVZ_DRbar(double);
   double calculate_MVWm_DRbar(double);

   double ThetaW() const;
   double v() const;
   double ThetaWp() const;

protected:
   std::size_t number_of_mass_iterations;
   unsigned pole_mass_loop_order;
   double diagonalization_precision;  ///< required precision in mass diagonalization
   bool calculate_sm_pole_masses; ///< switch to calculate the pole masses of the Standard Model particles
   bool force_output;             ///< switch to force output of pole masses
   CSE6SSM_physical physical; ///< contains the pole masses and mixings
   CSE6SSM_physical drbar;  ///< contains the running masses and mixings
   Two_loop_corrections two_loop_corrections; ///< used 2-loop corrections
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> problems;

   void copy_DRbar_masses_to_pole_masses();

private:

#ifdef ENABLE_THREADS
   struct Thread {
      typedef void(CSE6SSM_mass_eigenstates::*Memfun_t)();
      CSE6SSM_mass_eigenstates* model;
      Memfun_t fun;

      Thread(CSE6SSM_mass_eigenstates* model_, Memfun_t fun_)
         : model(model_), fun(fun_) {}
      void operator()() {
         try {
            (model->*fun)();
         } catch (...) {
            model->thread_exception = std::current_exception();
         }
      }
   };

   std::exception_ptr thread_exception;
   static std::mutex mtx_fortran; /// locks fortran functions
#endif

   int solve_ewsb_tree_level_via_soft_higgs_masses();

   // Passarino-Veltman loop functions
   double A0(double) const;
   double B0(double, double, double) const;
   double B1(double, double, double) const;
   double B00(double, double, double) const;
   double B22(double, double, double) const;
   double H0(double, double, double) const;
   double F0(double, double, double) const;
   double G0(double, double, double) const;

   // DR-bar masses
   double MVG;
   double MGlu;
   Eigen::Array<double,3,1> MFv;
   double MChaP;
   double MVP;
   double MVZ;
   double MVZp;
   Eigen::Array<double,6,1> MSd;
   Eigen::Array<double,3,1> MSv;
   Eigen::Array<double,6,1> MSu;
   Eigen::Array<double,6,1> MSe;
   Eigen::Array<double,6,1> MSDX;
   Eigen::Array<double,5,1> Mhh;
   Eigen::Array<double,5,1> MAh;
   Eigen::Array<double,2,1> MHpm;
   Eigen::Array<double,8,1> MChi;
   Eigen::Array<double,2,1> MCha;
   Eigen::Array<double,3,1> MFe;
   Eigen::Array<double,3,1> MFd;
   Eigen::Array<double,3,1> MFu;
   Eigen::Array<double,3,1> MFDX;
   Eigen::Array<double,7,1> MSHI0;
   Eigen::Array<double,4,1> MSHIPM;
   Eigen::Array<double,2,1> MChaI;
   Eigen::Array<double,7,1> MChiI;
   Eigen::Array<double,2,1> MSHp0;
   Eigen::Array<double,2,1> MSHpp;
   Eigen::Array<double,2,1> MChiP;
   double MVWm;

   // DR-bar mixing matrices
   Eigen::Matrix<double,6,6> ZD;
   Eigen::Matrix<double,3,3> ZV;
   Eigen::Matrix<double,6,6> ZU;
   Eigen::Matrix<double,6,6> ZE;
   Eigen::Matrix<double,6,6> ZDX;
   Eigen::Matrix<double,5,5> ZH;
   Eigen::Matrix<double,5,5> ZA;
   Eigen::Matrix<double,2,2> ZP;
   Eigen::Matrix<std::complex<double>,8,8> ZN;
   Eigen::Matrix<std::complex<double>,2,2> UM;
   Eigen::Matrix<std::complex<double>,2,2> UP;
   Eigen::Matrix<std::complex<double>,3,3> ZEL;
   Eigen::Matrix<std::complex<double>,3,3> ZER;
   Eigen::Matrix<std::complex<double>,3,3> ZDL;
   Eigen::Matrix<std::complex<double>,3,3> ZDR;
   Eigen::Matrix<std::complex<double>,3,3> ZUL;
   Eigen::Matrix<std::complex<double>,3,3> ZUR;
   Eigen::Matrix<std::complex<double>,3,3> ZDXL;
   Eigen::Matrix<std::complex<double>,3,3> ZDXR;
   Eigen::Matrix<double,7,7> UHI0;
   Eigen::Matrix<double,4,4> UHIPM;
   Eigen::Matrix<std::complex<double>,2,2> ZMI;
   Eigen::Matrix<std::complex<double>,2,2> ZPI;
   Eigen::Matrix<std::complex<double>,7,7> ZNI;
   Eigen::Matrix<double,2,2> UHp0;
   Eigen::Matrix<double,2,2> UHpp;
   Eigen::Matrix<std::complex<double>,2,2> ZNp;

   // phases
   std::complex<double> PhaseGlu;

};

std::ostream& operator<<(std::ostream&, const CSE6SSM_mass_eigenstates&);

} // namespace flexiblesusy

#endif
