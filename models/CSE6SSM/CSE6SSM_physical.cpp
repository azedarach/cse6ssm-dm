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

#include "CSE6SSM_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

CSE6SSM_physical::CSE6SSM_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MChaP(0), MVP(0),
       MVZ(0), MVZp(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array<
       double,3,1>::Zero()), MSu(Eigen::Array<double,6,1>::Zero()), MSe(
       Eigen::Array<double,6,1>::Zero()), MSDX(Eigen::Array<double,6,1>::Zero()),
       Mhh(Eigen::Array<double,5,1>::Zero()), MAh(Eigen::Array<double,5,1>::Zero()
       ), MHpm(Eigen::Array<double,2,1>::Zero()), MChi(Eigen::Array<double,8,1>
       ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,
       3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<
       double,3,1>::Zero()), MFDX(Eigen::Array<double,3,1>::Zero()), MSHI0(
       Eigen::Array<double,7,1>::Zero()), MSHIPM(Eigen::Array<double,4,1>::Zero())
       , MChaI(Eigen::Array<double,2,1>::Zero()), MChiI(Eigen::Array<double,7,1>
       ::Zero()), MSHp0(Eigen::Array<double,2,1>::Zero()), MSHpp(Eigen::Array<
       double,2,1>::Zero()), MChiP(Eigen::Array<double,2,1>::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZDX(Eigen::Matrix<double,6,6>::Zero()), ZH(Eigen::Matrix<double,5
      ,5>::Zero()), ZA(Eigen::Matrix<double,5,5>::Zero()), ZP(Eigen::Matrix<double
      ,2,2>::Zero()), ZN(Eigen::Matrix<std::complex<double>,8,8>::Zero()), UM(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDXL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDXR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), UHI0(Eigen::Matrix<double,7,7>::Zero()),
      UHIPM(Eigen::Matrix<double,4,4>::Zero()), ZMI(Eigen::Matrix<std::complex<
      double>,2,2>::Zero()), ZPI(Eigen::Matrix<std::complex<double>,2,2>::Zero()),
      ZNI(Eigen::Matrix<std::complex<double>,7,7>::Zero()), UHp0(Eigen::Matrix<
      double,2,2>::Zero()), UHpp(Eigen::Matrix<double,2,2>::Zero()), ZNp(
      Eigen::Matrix<std::complex<double>,2,2>::Zero())

{
}

void CSE6SSM_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MChaP = 0.;
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSDX = Eigen::Matrix<double,6,1>::Zero();
   ZDX = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,5,1>::Zero();
   ZH = Eigen::Matrix<double,5,5>::Zero();
   MAh = Eigen::Matrix<double,5,1>::Zero();
   ZA = Eigen::Matrix<double,5,5>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,8,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,8,8>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFDX = Eigen::Matrix<double,3,1>::Zero();
   ZDXL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDXR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MSHI0 = Eigen::Matrix<double,7,1>::Zero();
   UHI0 = Eigen::Matrix<double,7,7>::Zero();
   MSHIPM = Eigen::Matrix<double,4,1>::Zero();
   UHIPM = Eigen::Matrix<double,4,4>::Zero();
   MChaI = Eigen::Matrix<double,2,1>::Zero();
   ZMI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   ZPI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MChiI = Eigen::Matrix<double,7,1>::Zero();
   ZNI = Eigen::Matrix<std::complex<double>,7,7>::Zero();
   MSHp0 = Eigen::Matrix<double,2,1>::Zero();
   UHp0 = Eigen::Matrix<double,2,2>::Zero();
   MSHpp = Eigen::Matrix<double,2,1>::Zero();
   UHpp = Eigen::Matrix<double,2,2>::Zero();
   MChiP = Eigen::Matrix<double,2,1>::Zero();
   ZNp = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void CSE6SSM_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChiI), LOCALPHYSICAL(ZNI));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChiP), LOCALPHYSICAL(ZNp));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void CSE6SSM_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiI), LOCALPHYSICAL(ZNI));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiP), LOCALPHYSICAL(ZNp));

}

void CSE6SSM_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MChaP = " << MChaP << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MVZp = " << MVZp << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSDX = " << MSDX.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFDX = " << MFDX.transpose() << '\n';
   ostr << "MSHI0 = " << MSHI0.transpose() << '\n';
   ostr << "MSHIPM = " << MSHIPM.transpose() << '\n';
   ostr << "MChaI = " << MChaI.transpose() << '\n';
   ostr << "MChiI = " << MChiI.transpose() << '\n';
   ostr << "MSHp0 = " << MSHp0.transpose() << '\n';
   ostr << "MSHpp = " << MSHpp.transpose() << '\n';
   ostr << "MChiP = " << MChiP.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZDX = " << ZDX << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZDXL = " << ZDXL << '\n';
   ostr << "ZDXR = " << ZDXR << '\n';
   ostr << "UHI0 = " << UHI0 << '\n';
   ostr << "UHIPM = " << UHIPM << '\n';
   ostr << "ZMI = " << ZMI << '\n';
   ostr << "ZPI = " << ZPI << '\n';
   ostr << "ZNI = " << ZNI << '\n';
   ostr << "UHp0 = " << UHp0 << '\n';
   ostr << "UHpp = " << UHpp << '\n';
   ostr << "ZNp = " << ZNp << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const CSE6SSM_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
