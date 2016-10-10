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

// File generated at Wed 3 Jun 2015 23:47:48

#include "CSE6SSM_slha_io.hpp"
#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "CSE6SSM_two_scale_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALDRBAR(p) drbar.p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_DRBAR_PARAMETER(p) decltype(LOCALDRBAR(p)) p;
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p

using namespace softsusy;

namespace flexiblesusy {

char const * const CSE6SSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "NCharge", "Yu", "Yd", "Ye", "Te", "Td", "Tu", "HMIX", "ESIXHEYUK",
   "ESIXRUN", "ESIXGDYUK", "ESIXFUYUK", "ESIXFDYUK", "ESIXTHETRI", "ESIXTGDTRI"
   , "ESIXTFUTRI", "ESIXTFDTRI", "MSQ2", "MSE2", "MSL2", "MSU2", "MSD2",
   "MSOFT", "mX2", "mXBar2", "ESIXKAPPA", "ESIXTKAPPA", "ESIXLAMBDA",
   "ESIXTLAMBDA", "PHASES", "DRBARMASS", "DRBARUHNIMIX", "DRBARUHIPMIX",
   "DRBARUHNPMIX", "DRBARUHPPMIX", "DRBARUMIX", "DRBARVMIX", "DRBARNMAMIX",
   "DRBARDSQMIX", "DRBARESIXZDX", "DRBARESIXZXL", "DRBARESIXZXR",
   "DRBARSELMIX", "DRBARNMHMIX", "DRBARESIXZMI", "DRBARNMNMIX", "DRBARZNIMIX",
   "DRBARZNPMIX", "DRBARCHARGEMIX", "DRBARESIXZPI", "DRBARUSQMIX",
   "DRBARSNUMIX", "DRBARUELMIX", "DRBARUERMIX", "DRBARUDLMIX", "DRBARUDRMIX",
   "DRBARUULMIX", "DRBARUURMIX" }
;

CSE6SSM_slha_io::CSE6SSM_slha_io()
   : slha_io()
{
}

void CSE6SSM_slha_io::clear()
{
   slha_io.clear();
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CSE6SSM_slha_io::set_extpar(const CSE6SSM_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.sInput, "sInput");
   extpar << FORMAT_ELEMENT(72, input.QSInput, "QS");
   slha_io.set_block(extpar);

}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CSE6SSM_slha_io::set_extpar(const CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.sInput, "sInput");
   extpar << FORMAT_ELEMENT(72, input.QSInput, "QSInput");
   slha_io.set_block(extpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CSE6SSM_slha_io::set_minpar(const CSE6SSM_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.m0, "m0");
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(4, input.SignLambdax, "SignLambdax");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CSE6SSM_slha_io::set_minpar(const CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void CSE6SSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void CSE6SSM_slha_io::set_spinfo(const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_warning()) {
      std::ostringstream warnings;
      problems.print_warnings(warnings);
      spinfo << FORMAT_SPINFO(3, warnings.str());
   }

   if (problems.have_problem()) {
      std::ostringstream problems_str;
      problems.print_problems(problems_str);
      spinfo << FORMAT_SPINFO(4, problems_str.str());
   }

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the model DR-bar masses in the SLHA object
 *
 * @param drbar struct of DR-bar masses and mixings
 *
 * @param scale scale at which the DR-bar masses are evaluated
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void CSE6SSM_slha_io::set_drbar_mass(const CSE6SSM_physical& drbar,
                                     double scale,
                                     bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block DRBARMASS Q= " << FORMAT_SCALE(scale) << '\n'
      << FORMAT_MASS(1000021, LOCALDRBAR(MGlu), "DRbarGlu")
      << FORMAT_MASS(24, LOCALDRBAR(MVWm), "DRbarVWm")
      << FORMAT_MASS(1000091, LOCALDRBAR(MChaP), "DRbarChaP")
      << FORMAT_MASS(31, LOCALDRBAR(MVZp), "DRbarVZp")
      << FORMAT_MASS(1000092, LOCALDRBAR(MChiP(0)), "DRbarChiP(1)")
      << FORMAT_MASS(1000094, LOCALDRBAR(MChiP(1)), "DRbarChiP(2)")
      << FORMAT_MASS(1000024, LOCALDRBAR(MCha(0)), "DRbarCha(1)")
      << FORMAT_MASS(1000037, LOCALDRBAR(MCha(1)), "DRbarCha(2)")
      << FORMAT_MASS(37, LOCALDRBAR(MHpm(1)), "DRbarHpm(2)")
      << FORMAT_MASS(92, LOCALDRBAR(MSHp0(0)), "DRbarSHp0(1)")
      << FORMAT_MASS(94, LOCALDRBAR(MSHp0(1)), "DRbarSHp0(2)")
      << FORMAT_MASS(91, LOCALDRBAR(MSHpp(0)), "DRbarSHpp(1)")
      << FORMAT_MASS(93, LOCALDRBAR(MSHpp(1)), "DRbarSHpp(2)")
      << FORMAT_MASS(1000088, LOCALDRBAR(MChaI(0)), "DRbarChaI(1)")
      << FORMAT_MASS(1000089, LOCALDRBAR(MChaI(1)), "DRbarChaI(2)")
      << FORMAT_MASS(1000012, LOCALDRBAR(MSv(0)), "DRbarSv(1)")
      << FORMAT_MASS(1000014, LOCALDRBAR(MSv(1)), "DRbarSv(2)")
      << FORMAT_MASS(1000016, LOCALDRBAR(MSv(2)), "DRbarSv(3)")
      << FORMAT_MASS(51, LOCALDRBAR(MFDX(0)), "DRbarFDX(1)")
      << FORMAT_MASS(52, LOCALDRBAR(MFDX(1)), "DRbarFDX(2)")
      << FORMAT_MASS(53, LOCALDRBAR(MFDX(2)), "DRbarFDX(3)")
      << FORMAT_MASS(81, LOCALDRBAR(MSHIPM(0)), "DRbarSHIPM(1)")
      << FORMAT_MASS(85, LOCALDRBAR(MSHIPM(1)), "DRbarSHIPM(2)")
      << FORMAT_MASS(83, LOCALDRBAR(MSHIPM(2)), "DRbarSHIPM(3)")
      << FORMAT_MASS(87, LOCALDRBAR(MSHIPM(3)), "DRbarSHIPM(4)")
      << FORMAT_MASS(25, LOCALDRBAR(Mhh(0)), "DRbarhh(1)")
      << FORMAT_MASS(35, LOCALDRBAR(Mhh(1)), "DRbarhh(2)")
      << FORMAT_MASS(45, LOCALDRBAR(Mhh(2)), "DRbarhh(3)")
      << FORMAT_MASS(55, LOCALDRBAR(Mhh(3)), "DRbarhh(4)")
      << FORMAT_MASS(65, LOCALDRBAR(Mhh(4)), "DRbarhh(5)")
      << FORMAT_MASS(91191138, LOCALDRBAR(MAh(2)), "DRbarAh(3)")
      << FORMAT_MASS(36, LOCALDRBAR(MAh(3)), "DRbarAh(4)")
      << FORMAT_MASS(91191137, LOCALDRBAR(MAh(4)), "DRbarAh(5)")
      << FORMAT_MASS(1000001, LOCALDRBAR(MSd(0)), "DRbarSd(1)")
      << FORMAT_MASS(1000003, LOCALDRBAR(MSd(1)), "DRbarSd(2)")
      << FORMAT_MASS(1000005, LOCALDRBAR(MSd(2)), "DRbarSd(3)")
      << FORMAT_MASS(2000001, LOCALDRBAR(MSd(3)), "DRbarSd(4)")
      << FORMAT_MASS(2000003, LOCALDRBAR(MSd(4)), "DRbarSd(5)")
      << FORMAT_MASS(2000005, LOCALDRBAR(MSd(5)), "DRbarSd(6)")
      << FORMAT_MASS(1000011, LOCALDRBAR(MSe(0)), "DRbarSe(1)")
      << FORMAT_MASS(1000013, LOCALDRBAR(MSe(1)), "DRbarSe(2)")
      << FORMAT_MASS(1000015, LOCALDRBAR(MSe(2)), "DRbarSe(3)")
      << FORMAT_MASS(2000011, LOCALDRBAR(MSe(3)), "DRbarSe(4)")
      << FORMAT_MASS(2000013, LOCALDRBAR(MSe(4)), "DRbarSe(5)")
      << FORMAT_MASS(2000015, LOCALDRBAR(MSe(5)), "DRbarSe(6)")
      << FORMAT_MASS(1000002, LOCALDRBAR(MSu(0)), "DRbarSu(1)")
      << FORMAT_MASS(1000004, LOCALDRBAR(MSu(1)), "DRbarSu(2)")
      << FORMAT_MASS(1000006, LOCALDRBAR(MSu(2)), "DRbarSu(3)")
      << FORMAT_MASS(2000002, LOCALDRBAR(MSu(3)), "DRbarSu(4)")
      << FORMAT_MASS(2000004, LOCALDRBAR(MSu(4)), "DRbarSu(5)")
      << FORMAT_MASS(2000006, LOCALDRBAR(MSu(5)), "DRbarSu(6)")
      << FORMAT_MASS(1000051, LOCALDRBAR(MSDX(0)), "DRbarSDX(1)")
      << FORMAT_MASS(2000051, LOCALDRBAR(MSDX(1)), "DRbarSDX(2)")
      << FORMAT_MASS(1000052, LOCALDRBAR(MSDX(2)), "DRbarSDX(3)")
      << FORMAT_MASS(2000052, LOCALDRBAR(MSDX(3)), "DRbarSDX(4)")
      << FORMAT_MASS(1000053, LOCALDRBAR(MSDX(4)), "DRbarSDX(5)")
      << FORMAT_MASS(2000053, LOCALDRBAR(MSDX(5)), "DRbarSDX(6)")
      << FORMAT_MASS(1000081, LOCALDRBAR(MChiI(0)), "DRbarChiI(1)")
      << FORMAT_MASS(1000082, LOCALDRBAR(MChiI(1)), "DRbarChiI(2)")
      << FORMAT_MASS(1000083, LOCALDRBAR(MChiI(2)), "DRbarChiI(3)")
      << FORMAT_MASS(1000084, LOCALDRBAR(MChiI(3)), "DRbarChiI(4)")
      << FORMAT_MASS(1000085, LOCALDRBAR(MChiI(4)), "DRbarChiI(5)")
      << FORMAT_MASS(1000086, LOCALDRBAR(MChiI(5)), "DRbarChiI(6)")
      << FORMAT_MASS(1000087, LOCALDRBAR(MChiI(6)), "DRbarChiI(7)")
      << FORMAT_MASS(82, LOCALDRBAR(MSHI0(0)), "DRbarSHI0(1)")
      << FORMAT_MASS(86, LOCALDRBAR(MSHI0(1)), "DRbarSHI0(2)")
      << FORMAT_MASS(84, LOCALDRBAR(MSHI0(2)), "DRbarSHI0(3)")
      << FORMAT_MASS(88, LOCALDRBAR(MSHI0(3)), "DRbarSHI0(4)")
      << FORMAT_MASS(9994453, LOCALDRBAR(MSHI0(4)), "DRbarSHI0(5)")
      << FORMAT_MASS(9994454, LOCALDRBAR(MSHI0(5)), "DRbarSHI0(6)")
      << FORMAT_MASS(9994455, LOCALDRBAR(MSHI0(6)), "DRbarSHI0(7)")
      << FORMAT_MASS(1000022, LOCALDRBAR(MChi(0)), "DRbarChi(1)")
      << FORMAT_MASS(1000023, LOCALDRBAR(MChi(1)), "DRbarChi(2)")
      << FORMAT_MASS(1000025, LOCALDRBAR(MChi(2)), "DRbarChi(3)")
      << FORMAT_MASS(1000035, LOCALDRBAR(MChi(3)), "DRbarChi(4)")
      << FORMAT_MASS(1000045, LOCALDRBAR(MChi(4)), "DRbarChi(5)")
      << FORMAT_MASS(1000055, LOCALDRBAR(MChi(5)), "DRbarChi(6)")
      << FORMAT_MASS(1000065, LOCALDRBAR(MChi(6)), "DRbarChi(7)")
      << FORMAT_MASS(1000075, LOCALDRBAR(MChi(7)), "DRbarChi(8)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALDRBAR(MVG), "DRbarVG")
         << FORMAT_MASS(12, LOCALDRBAR(MFv(0)), "DRbarFv(1)")
         << FORMAT_MASS(14, LOCALDRBAR(MFv(1)), "DRbarFv(2)")
         << FORMAT_MASS(16, LOCALDRBAR(MFv(2)), "DRbarFv(3)")
         << FORMAT_MASS(22, LOCALDRBAR(MVP), "DRbarVP")
         << FORMAT_MASS(23, LOCALDRBAR(MVZ), "DRbarVZ")
         << FORMAT_MASS(11, LOCALDRBAR(MFe(0)), "DRbarFe(1)")
         << FORMAT_MASS(13, LOCALDRBAR(MFe(1)), "DRbarFe(2)")
         << FORMAT_MASS(15, LOCALDRBAR(MFe(2)), "DRbarFe(3)")
         << FORMAT_MASS(1, LOCALDRBAR(MFd(0)), "DRbarFd(1)")
         << FORMAT_MASS(3, LOCALDRBAR(MFd(1)), "DRbarFd(2)")
         << FORMAT_MASS(5, LOCALDRBAR(MFd(2)), "DRbarFd(3)")
         << FORMAT_MASS(2, LOCALDRBAR(MFu(0)), "DRbarFu(1)")
         << FORMAT_MASS(4, LOCALDRBAR(MFu(1)), "DRbarFu(2)")
         << FORMAT_MASS(6, LOCALDRBAR(MFu(2)), "DRbarFu(3)")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the DR-bar mixing matrices in the SLHA object.
 *
 * @param physical struct of DR-bar parameters
 *
 * @param scale scale at which the DR-bar mixings are evaluated
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void CSE6SSM_slha_io::set_drbar_mixing_matrices(const CSE6SSM_physical& drbar,
                                                double scale,
                                                bool write_sm_mixing_matrics)
{
   slha_io.set_block("DRBARUHNIMIX", LOCALDRBAR(UHI0), "DRbarUHI0", scale);
   slha_io.set_block("DRBARUHIPMIX", LOCALDRBAR(UHIPM), "DRbarUHIPM", scale);
   slha_io.set_block("DRBARUHNPMIX", LOCALDRBAR(UHp0), "DRbarUHp0", scale);
   slha_io.set_block("DRBARUHPPMIX", LOCALDRBAR(UHpp), "DRbarUHpp", scale);
   slha_io.set_block("DRBARUMIX", LOCALDRBAR(UM), "DRbarUM", scale);
   slha_io.set_block("DRBARVMIX", LOCALDRBAR(UP), "DRbarUP", scale);
   slha_io.set_block("DRBARNMAMIX", LOCALDRBAR(ZA), "DRbarZA", scale);
   slha_io.set_block("DRBARDSQMIX", LOCALDRBAR(ZD), "DRbarZD", scale);
   slha_io.set_block("DRBARESIXZDX", LOCALDRBAR(ZDX), "DRbarZDX", scale);
   slha_io.set_block("DRBARESIXZXL", LOCALDRBAR(ZDXL), "DRbarZDXL", scale);
   slha_io.set_block("DRBARESIXZXR", LOCALDRBAR(ZDXR), "DRbarZDXR", scale);
   slha_io.set_block("DRBARSELMIX", LOCALDRBAR(ZE), "DRbarZE", scale);
   slha_io.set_block("DRBARNMHMIX", LOCALDRBAR(ZH), "DRbarZH", scale);
   slha_io.set_block("DRBARESIXZMI", LOCALDRBAR(ZMI), "DRbarZMI", scale);
   slha_io.set_block("DRBARNMNMIX", LOCALDRBAR(ZN), "DRbarZN", scale);
   slha_io.set_block("DRBARZNIMIX", LOCALDRBAR(ZNI), "DRbarZNI", scale);
   slha_io.set_block("DRBARZNPMIX", LOCALDRBAR(ZNp), "DRbarZNp", scale);
   slha_io.set_block("DRBARCHARGEMIX", LOCALDRBAR(ZP), "DRbarZP", scale);
   slha_io.set_block("DRBARESIXZPI", LOCALDRBAR(ZPI), "DRbarZPI", scale);
   slha_io.set_block("DRBARUSQMIX", LOCALDRBAR(ZU), "DRbarZU", scale);
   slha_io.set_block("DRBARSNUMIX", LOCALDRBAR(ZV), "DRbarZV", scale);

   if (write_sm_mixing_matrics) {
      slha_io.set_block("DRBARUELMIX", LOCALDRBAR(ZEL), "DRbarZEL", scale);
      slha_io.set_block("DRBARUERMIX", LOCALDRBAR(ZER), "DRbarZER", scale);
      slha_io.set_block("DRBARUDLMIX", LOCALDRBAR(ZDL), "DRbarZDL", scale);
      slha_io.set_block("DRBARUDRMIX", LOCALDRBAR(ZDR), "DRbarZDR", scale);
      slha_io.set_block("DRBARUULMIX", LOCALDRBAR(ZUL), "DRbarZUL", scale);
      slha_io.set_block("DRBARUURMIX", LOCALDRBAR(ZUR), "DRbarZUR", scale);
   }

}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void CSE6SSM_slha_io::set_mass(const CSE6SSM_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(1000091, LOCALPHYSICAL(MChaP), "ChaP")
      << FORMAT_MASS(31, LOCALPHYSICAL(MVZp), "VZp")
      << FORMAT_MASS(1000092, LOCALPHYSICAL(MChiP(0)), "ChiP(1)")
      << FORMAT_MASS(1000094, LOCALPHYSICAL(MChiP(1)), "ChiP(2)")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha(2)")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHpm(1)), "Hpm(2)")
      << FORMAT_MASS(92, LOCALPHYSICAL(MSHp0(0)), "SHp0(1)")
      << FORMAT_MASS(94, LOCALPHYSICAL(MSHp0(1)), "SHp0(2)")
      << FORMAT_MASS(91, LOCALPHYSICAL(MSHpp(0)), "SHpp(1)")
      << FORMAT_MASS(93, LOCALPHYSICAL(MSHpp(1)), "SHpp(2)")
      << FORMAT_MASS(1000088, LOCALPHYSICAL(MChaI(0)), "ChaI(1)")
      << FORMAT_MASS(1000089, LOCALPHYSICAL(MChaI(1)), "ChaI(2)")
      << FORMAT_MASS(1000012, LOCALPHYSICAL(MSv(0)), "Sv(1)")
      << FORMAT_MASS(1000014, LOCALPHYSICAL(MSv(1)), "Sv(2)")
      << FORMAT_MASS(1000016, LOCALPHYSICAL(MSv(2)), "Sv(3)")
      << FORMAT_MASS(51, LOCALPHYSICAL(MFDX(0)), "FDX(1)")
      << FORMAT_MASS(52, LOCALPHYSICAL(MFDX(1)), "FDX(2)")
      << FORMAT_MASS(53, LOCALPHYSICAL(MFDX(2)), "FDX(3)")
      << FORMAT_MASS(81, LOCALPHYSICAL(MSHIPM(0)), "SHIPM(1)")
      << FORMAT_MASS(85, LOCALPHYSICAL(MSHIPM(1)), "SHIPM(2)")
      << FORMAT_MASS(83, LOCALPHYSICAL(MSHIPM(2)), "SHIPM(3)")
      << FORMAT_MASS(87, LOCALPHYSICAL(MSHIPM(3)), "SHIPM(4)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(45, LOCALPHYSICAL(Mhh(2)), "hh(3)")
      << FORMAT_MASS(55, LOCALPHYSICAL(Mhh(3)), "hh(4)")
      << FORMAT_MASS(65, LOCALPHYSICAL(Mhh(4)), "hh(5)")
      << FORMAT_MASS(91191138, LOCALPHYSICAL(MAh(2)), "Ah(3)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(3)), "Ah(4)")
      << FORMAT_MASS(91191137, LOCALPHYSICAL(MAh(4)), "Ah(5)")
      << FORMAT_MASS(1000001, LOCALPHYSICAL(MSd(0)), "Sd(1)")
      << FORMAT_MASS(1000003, LOCALPHYSICAL(MSd(1)), "Sd(2)")
      << FORMAT_MASS(1000005, LOCALPHYSICAL(MSd(2)), "Sd(3)")
      << FORMAT_MASS(2000001, LOCALPHYSICAL(MSd(3)), "Sd(4)")
      << FORMAT_MASS(2000003, LOCALPHYSICAL(MSd(4)), "Sd(5)")
      << FORMAT_MASS(2000005, LOCALPHYSICAL(MSd(5)), "Sd(6)")
      << FORMAT_MASS(1000011, LOCALPHYSICAL(MSe(0)), "Se(1)")
      << FORMAT_MASS(1000013, LOCALPHYSICAL(MSe(1)), "Se(2)")
      << FORMAT_MASS(1000015, LOCALPHYSICAL(MSe(2)), "Se(3)")
      << FORMAT_MASS(2000011, LOCALPHYSICAL(MSe(3)), "Se(4)")
      << FORMAT_MASS(2000013, LOCALPHYSICAL(MSe(4)), "Se(5)")
      << FORMAT_MASS(2000015, LOCALPHYSICAL(MSe(5)), "Se(6)")
      << FORMAT_MASS(1000002, LOCALPHYSICAL(MSu(0)), "Su(1)")
      << FORMAT_MASS(1000004, LOCALPHYSICAL(MSu(1)), "Su(2)")
      << FORMAT_MASS(1000006, LOCALPHYSICAL(MSu(2)), "Su(3)")
      << FORMAT_MASS(2000002, LOCALPHYSICAL(MSu(3)), "Su(4)")
      << FORMAT_MASS(2000004, LOCALPHYSICAL(MSu(4)), "Su(5)")
      << FORMAT_MASS(2000006, LOCALPHYSICAL(MSu(5)), "Su(6)")
      << FORMAT_MASS(1000051, LOCALPHYSICAL(MSDX(0)), "SDX(1)")
      << FORMAT_MASS(2000051, LOCALPHYSICAL(MSDX(1)), "SDX(2)")
      << FORMAT_MASS(1000052, LOCALPHYSICAL(MSDX(2)), "SDX(3)")
      << FORMAT_MASS(2000052, LOCALPHYSICAL(MSDX(3)), "SDX(4)")
      << FORMAT_MASS(1000053, LOCALPHYSICAL(MSDX(4)), "SDX(5)")
      << FORMAT_MASS(2000053, LOCALPHYSICAL(MSDX(5)), "SDX(6)")
      << FORMAT_MASS(1000081, LOCALPHYSICAL(MChiI(0)), "ChiI(1)")
      << FORMAT_MASS(1000082, LOCALPHYSICAL(MChiI(1)), "ChiI(2)")
      << FORMAT_MASS(1000083, LOCALPHYSICAL(MChiI(2)), "ChiI(3)")
      << FORMAT_MASS(1000084, LOCALPHYSICAL(MChiI(3)), "ChiI(4)")
      << FORMAT_MASS(1000085, LOCALPHYSICAL(MChiI(4)), "ChiI(5)")
      << FORMAT_MASS(1000086, LOCALPHYSICAL(MChiI(5)), "ChiI(6)")
      << FORMAT_MASS(1000087, LOCALPHYSICAL(MChiI(6)), "ChiI(7)")
      << FORMAT_MASS(82, LOCALPHYSICAL(MSHI0(0)), "SHI0(1)")
      << FORMAT_MASS(86, LOCALPHYSICAL(MSHI0(1)), "SHI0(2)")
      << FORMAT_MASS(84, LOCALPHYSICAL(MSHI0(2)), "SHI0(3)")
      << FORMAT_MASS(88, LOCALPHYSICAL(MSHI0(3)), "SHI0(4)")
      << FORMAT_MASS(9994453, LOCALPHYSICAL(MSHI0(4)), "SHI0(5)")
      << FORMAT_MASS(9994454, LOCALPHYSICAL(MSHI0(5)), "SHI0(6)")
      << FORMAT_MASS(9994455, LOCALPHYSICAL(MSHI0(6)), "SHI0(7)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
      << FORMAT_MASS(1000045, LOCALPHYSICAL(MChi(4)), "Chi(5)")
      << FORMAT_MASS(1000055, LOCALPHYSICAL(MChi(5)), "Chi(6)")
      << FORMAT_MASS(1000065, LOCALPHYSICAL(MChi(6)), "Chi(7)")
      << FORMAT_MASS(1000075, LOCALPHYSICAL(MChi(7)), "Chi(8)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void CSE6SSM_slha_io::set_mixing_matrices(const CSE6SSM_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UHNIMIX", LOCALPHYSICAL(UHI0), "UHI0");
   slha_io.set_block("UHIPMIX", LOCALPHYSICAL(UHIPM), "UHIPM");
   slha_io.set_block("UHNPMIX", LOCALPHYSICAL(UHp0), "UHp0");
   slha_io.set_block("UHPPMIX", LOCALPHYSICAL(UHpp), "UHpp");
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("NMAMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", LOCALPHYSICAL(ZD), "ZD");
   slha_io.set_block("ESIXZDX", LOCALPHYSICAL(ZDX), "ZDX");
   slha_io.set_block("ESIXZXL", LOCALPHYSICAL(ZDXL), "ZDXL");
   slha_io.set_block("ESIXZXR", LOCALPHYSICAL(ZDXR), "ZDXR");
   slha_io.set_block("SELMIX", LOCALPHYSICAL(ZE), "ZE");
   slha_io.set_block("NMHMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("ESIXZMI", LOCALPHYSICAL(ZMI), "ZMI");
   slha_io.set_block("NMNMIX", LOCALPHYSICAL(ZN), "ZN");
   slha_io.set_block("ZNIMIX", LOCALPHYSICAL(ZNI), "ZNI");
   slha_io.set_block("ZNPMIX", LOCALPHYSICAL(ZNp), "ZNp");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");
   slha_io.set_block("ESIXZPI", LOCALPHYSICAL(ZPI), "ZPI");
   slha_io.set_block("USQMIX", LOCALPHYSICAL(ZU), "ZU");
   slha_io.set_block("SNUMIX", LOCALPHYSICAL(ZV), "ZV");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UELMIX", LOCALPHYSICAL(ZEL), "ZEL");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(ZER), "ZER");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(ZDL), "ZDL");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(ZDR), "ZDR");
      slha_io.set_block("UULMIX", LOCALPHYSICAL(ZUL), "ZUL");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(ZUR), "ZUR");
   }

}

void CSE6SSM_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void CSE6SSM_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Write SLHA object to file.
 *
 * @param file_name file name
 */
void CSE6SSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double CSE6SSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void CSE6SSM_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void CSE6SSM_slha_io::fill(CSE6SSM_input_parameters<Two_scale>& input) const
{
   void (*fill_two_scale_minpar_tuple) (CSE6SSM_input_parameters<Two_scale>&,
                                        int, double) = &CSE6SSM_slha_io::fill_minpar_tuple;
   void (*fill_two_scale_extpar_tuple) (CSE6SSM_input_parameters<Two_scale>&,
                                        int, double) = &CSE6SSM_slha_io::fill_extpar_tuple;

   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_two_scale_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_two_scale_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   input.MuPhiInput = slha_io.read_entry("HMIXIN", 31);
   input.KappaPrInput = slha_io.read_entry("HMIXIN", 32);
   input.SigmaxInput = slha_io.read_entry("HMIXIN", 33);
   slha_io.read_block("ESIXHEYUKIN", input.hEInput);
   input.SigmaLInput = slha_io.read_entry("ESIXRUNIN", 42);
   slha_io.read_block("ESIXGDYUKIN", input.gDInput);
   slha_io.read_block("ESIXFUYUKIN", input.fuInput);
   slha_io.read_block("ESIXFDYUKIN", input.fdInput);
   input.BMuPhiInput = slha_io.read_entry("ESIXRUNIN", 30);
   slha_io.read_block("ESIXKAPPAIN", input.KappaInput);
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.BMuPrInput = slha_io.read_entry("ESIXRUNIN", 101);

}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void CSE6SSM_slha_io::fill(CSE6SSM_semianalytic_input_parameters<Two_scale>& input) const
{
   void (*fill_semi_two_scale_minpar_tuple) (CSE6SSM_semianalytic_input_parameters<Two_scale>&,
                                        int, double) = &CSE6SSM_slha_io::fill_minpar_tuple;
   void (*fill_semi_two_scale_extpar_tuple) (CSE6SSM_semianalytic_input_parameters<Two_scale>&,
                                        int, double) = &CSE6SSM_slha_io::fill_extpar_tuple;

   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_semi_two_scale_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_semi_two_scale_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   input.LambdaxInput = slha_io.read_entry("ESIXRUNIN", 1); 
   input.MuPhiInput = slha_io.read_entry("HMIXIN", 31);
   input.KappaPrInput = slha_io.read_entry("HMIXIN", 32);
   input.SigmaxInput = slha_io.read_entry("HMIXIN", 33);
   slha_io.read_block("ESIXHEYUKIN", input.hEInput);
   input.SigmaLInput = slha_io.read_entry("ESIXRUNIN", 42);
   slha_io.read_block("ESIXGDYUKIN", input.gDInput);
   slha_io.read_block("ESIXFUYUKIN", input.fuInput);
   slha_io.read_block("ESIXFDYUKIN", input.fdInput);
   input.BMuPhiInput = slha_io.read_entry("ESIXRUNIN", 30);
   slha_io.read_block("ESIXKAPPAIN", input.KappaInput);
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.BMuPrInput = slha_io.read_entry("ESIXRUNIN", 101);

}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
void CSE6SSM_slha_io::fill_drbar_parameters(CSE6SSM_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   model.set_g1p(slha_io.read_entry("gauge", 4));
   model.set_QS(slha_io.read_entry("NCharge", 1));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   {
      DEFINE_PARAMETER(TYe);
      slha_io.read_block("Te", TYe);
      model.set_TYe(TYe);
   }
   {
      DEFINE_PARAMETER(TYd);
      slha_io.read_block("Td", TYd);
      model.set_TYd(TYd);
   }
   {
      DEFINE_PARAMETER(TYu);
      slha_io.read_block("Tu", TYu);
      model.set_TYu(TYu);
   }
   model.set_XiF(slha_io.read_entry("HMIX", 30));
   model.set_MuPhi(slha_io.read_entry("HMIX", 31));
   model.set_KappaPr(slha_io.read_entry("HMIX", 32));
   model.set_Sigmax(slha_io.read_entry("HMIX", 33));
   {
      DEFINE_PARAMETER(hE);
      slha_io.read_block("ESIXHEYUK", hE);
      model.set_hE(hE);
   }
   model.set_SigmaL(slha_io.read_entry("ESIXRUN", 42));
   {
      DEFINE_PARAMETER(gD);
      slha_io.read_block("ESIXGDYUK", gD);
      model.set_gD(gD);
   }
   {
      DEFINE_PARAMETER(fu);
      slha_io.read_block("ESIXFUYUK", fu);
      model.set_fu(fu);
   }
   {
      DEFINE_PARAMETER(fd);
      slha_io.read_block("ESIXFDYUK", fd);
      model.set_fd(fd);
   }
   model.set_TKappaPr(slha_io.read_entry("ESIXRUN", 28));
   model.set_TSigmax(slha_io.read_entry("ESIXRUN", 29));
   {
      DEFINE_PARAMETER(ThE);
      slha_io.read_block("ESIXTHETRI", ThE);
      model.set_ThE(ThE);
   }
   model.set_TSigmaL(slha_io.read_entry("ESIXRUN", 43));
   {
      DEFINE_PARAMETER(TgD);
      slha_io.read_block("ESIXTGDTRI", TgD);
      model.set_TgD(TgD);
   }
   {
      DEFINE_PARAMETER(Tfu);
      slha_io.read_block("ESIXTFUTRI", Tfu);
      model.set_Tfu(Tfu);
   }
   {
      DEFINE_PARAMETER(Tfd);
      slha_io.read_block("ESIXTFDTRI", Tfd);
      model.set_Tfd(Tfd);
   }
   model.set_BMuPhi(slha_io.read_entry("ESIXRUN", 30));
   model.set_LXiF(slha_io.read_entry("HMIX", 34));
   {
      DEFINE_PARAMETER(mq2);
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      DEFINE_PARAMETER(me2);
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      DEFINE_PARAMETER(ml2);
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      DEFINE_PARAMETER(mu2);
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      DEFINE_PARAMETER(md2);
      slha_io.read_block("MSD2", md2);
      model.set_md2(md2);
   }
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   {
      DEFINE_PARAMETER(mDx2);
      slha_io.read_block("mX2", mDx2);
      model.set_mDx2(mDx2);
   }
   {
      DEFINE_PARAMETER(mDxbar2);
      slha_io.read_block("mXBar2", mDxbar2);
      model.set_mDxbar2(mDxbar2);
   }
   model.set_ms2(slha_io.read_entry("MSOFT", 23));
   model.set_msbar2(slha_io.read_entry("MSOFT", 24));
   model.set_mphi2(slha_io.read_entry("MSOFT", 25));
   model.set_mHp2(slha_io.read_entry("MSOFT", 26));
   model.set_mHpbar2(slha_io.read_entry("MSOFT", 27));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_MassBp(slha_io.read_entry("MSOFT", 4));
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));
   model.set_vs(slha_io.read_entry("ESIXRUN", 11));
   model.set_vsb(slha_io.read_entry("ESIXRUN", 12));
   model.set_vphi(slha_io.read_entry("ESIXRUN", 13));
   model.set_QS(slha_io.read_entry("EXTPAR", 72));
   {
      DEFINE_PARAMETER(Kappa);
      slha_io.read_block("ESIXKAPPA", Kappa);
      model.set_Kappa(Kappa);
   }
   {
      DEFINE_PARAMETER(TKappa);
      slha_io.read_block("ESIXTKAPPA", TKappa);
      model.set_TKappa(TKappa);
   }
   model.set_Lambdax(slha_io.read_entry("ESIXRUN", 1));
   model.set_TLambdax(slha_io.read_entry("ESIXRUN", 2));
   {
      DEFINE_PARAMETER(Lambda12);
      slha_io.read_block("ESIXLAMBDA", Lambda12);
      model.set_Lambda12(Lambda12);
   }
   {
      DEFINE_PARAMETER(TLambda12);
      slha_io.read_block("ESIXTLAMBDA", TLambda12);
      model.set_TLambda12(TLambda12);
   }
   model.set_MuPr(slha_io.read_entry("ESIXRUN", 0));
   model.set_BMuPr(slha_io.read_entry("ESIXRUN", 101));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 */
void CSE6SSM_slha_io::fill(CSE6SSM_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   CSE6SSM_physical drbar_hk;
   fill_drbar(drbar_hk);
   drbar_hk.convert_to_hk();
   model.get_drbar_masses() = drbar_hk;

   CSE6SSM_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void CSE6SSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&CSE6SSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void CSE6SSM_slha_io::fill_minpar_tuple(CSE6SSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.m0 = value; break;
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 4: input.SignLambdax = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CSE6SSM_slha_io::fill_minpar_tuple(CSE6SSM_semianalytic_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CSE6SSM_slha_io::fill_extpar_tuple(CSE6SSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 65: input.sInput = value; break;
   case 72: input.QSInput = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CSE6SSM_slha_io::fill_extpar_tuple(CSE6SSM_semianalytic_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 65: input.sInput = value; break;
   case 72: input.QSInput = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CSE6SSM_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
                                                  int key, double value)
{
   if (0 <= key && key < static_cast<int>(Spectrum_generator_settings::NUMBER_OF_OPTIONS)) {
      settings.set((Spectrum_generator_settings::Settings)key, value);
   } else {
      WARNING("Unrecognized key in block FlexibleSUSY: " << key);
   }
}

/**
 * Reads DR-bar masses and mixing matrices from a SLHA output file.
 */
void CSE6SSM_slha_io::fill_drbar(CSE6SSM_physical& drbar) const
{
   {
      DEFINE_DRBAR_PARAMETER(ZD);
      slha_io.read_block("DRBARDSQMIX", ZD);
      LOCALDRBAR(ZD) = ZD;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZU);
      slha_io.read_block("DRBARUSQMIX", ZU);
      LOCALDRBAR(ZU) = ZU;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZE);
      slha_io.read_block("DRBARSELMIX", ZE);
      LOCALDRBAR(ZE) = ZE;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZV);
      slha_io.read_block("DRBARSNUMIX", ZV);
      LOCALDRBAR(ZV) = ZV;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZDX);
      slha_io.read_block("DRBARESIXZDX", ZDX);
      LOCALDRBAR(ZDX) = ZDX;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZH);
      slha_io.read_block("DRBARNMHMIX", ZH);
      LOCALDRBAR(ZH) = ZH;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZA);
      slha_io.read_block("DRBARNMAMIX", ZA);
      LOCALDRBAR(ZA) = ZA;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZP);
      slha_io.read_block("DRBARCHARGEMIX", ZP);
      LOCALDRBAR(ZP) = ZP;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZN);
      slha_io.read_block("DRBARNMNMIX", ZN);
      LOCALDRBAR(ZN) = ZN;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZNp);
      slha_io.read_block("DRBARZNPMIX", ZNp);
      LOCALDRBAR(ZNp) = ZNp;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZNI);
      slha_io.read_block("DRBARMIX", ZNI);
      LOCALDRBAR(ZNI) = ZNI;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZMI);
      slha_io.read_block("DRBARESIXZMI", ZMI);
      LOCALDRBAR(ZMI) = ZMI;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZPI);
      slha_io.read_block("DRBARESIXZPI", ZPI);
      LOCALDRBAR(ZPI) = ZPI;
   }
   {
      DEFINE_DRBAR_PARAMETER(UP);
      slha_io.read_block("DRBARVMIX", UP);
      LOCALDRBAR(UP) = UP;
   }
   {
      DEFINE_DRBAR_PARAMETER(UM);
      slha_io.read_block("DRBARUMIX", UM);
      LOCALDRBAR(UM) = UM;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZEL);
      slha_io.read_block("DRBARUELMIX", ZEL);
      LOCALDRBAR(ZEL) = ZEL;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZER);
      slha_io.read_block("DRBARUERMIX", ZER);
      LOCALDRBAR(ZER) = ZER;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZDL);
      slha_io.read_block("DRBARUDLMIX", ZDL);
      LOCALDRBAR(ZDL) = ZDL;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZDR);
      slha_io.read_block("DRBARUDRMIX", ZDR);
      LOCALDRBAR(ZDR) = ZDR;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZUL);
      slha_io.read_block("DRBARUULMIX", ZUL);
      LOCALDRBAR(ZUL) = ZUL;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZUR);
      slha_io.read_block("DRBARUURMIX", ZUR);
      LOCALDRBAR(ZUR) = ZUR;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZDXL);
      slha_io.read_block("DRBARESIXZXL", ZDXL);
      LOCALDRBAR(ZDXL) = ZDXL;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZDXR);
      slha_io.read_block("DRBARESIXZXR", ZDXR);
      LOCALDRBAR(ZDXR) = ZDXR;
   }
   {
      DEFINE_DRBAR_PARAMETER(UHp0);
      slha_io.read_block("DRBARUHNPMIX", UHp0);
      LOCALDRBAR(UHp0) = UHp0;
   }
   {
      DEFINE_DRBAR_PARAMETER(UHpp);
      slha_io.read_block("DRBARUHPPMIX", UHpp);
      LOCALDRBAR(UHpp) = UHpp;
   }
   {
      DEFINE_DRBAR_PARAMETER(UHI0);
      slha_io.read_block("DRBARUHNIMIX", UHI0);
      LOCALDRBAR(UHI0) = UHI0;
   }
   {
      DEFINE_DRBAR_PARAMETER(UHIPM);
      slha_io.read_block("DRBARUHPPMIX", UHIPM);
      LOCALDRBAR(UHIPM) = UHIPM;
   }

   LOCALDRBAR(MVG) = slha_io.read_entry("DRBARMASS", 21);
   LOCALDRBAR(MGlu) = slha_io.read_entry("DRBARMASS", 1000021);
   LOCALDRBAR(MFv)(0) = slha_io.read_entry("DRBARMASS", 12);
   LOCALDRBAR(MFv)(1) = slha_io.read_entry("DRBARMASS", 14);
   LOCALDRBAR(MFv)(2) = slha_io.read_entry("DRBARMASS", 16);
   LOCALDRBAR(MChaP) = slha_io.read_entry("DRBARMASS", 1000091);
   LOCALDRBAR(MVP) = slha_io.read_entry("DRBARMASS", 22);
   LOCALDRBAR(MVZ) = slha_io.read_entry("DRBARMASS", 23);
   LOCALDRBAR(MVZp) = slha_io.read_entry("DRBARMASS", 31);
   LOCALDRBAR(MSd)(0) = slha_io.read_entry("DRBARMASS", 1000001);
   LOCALDRBAR(MSd)(1) = slha_io.read_entry("DRBARMASS", 1000003);
   LOCALDRBAR(MSd)(2) = slha_io.read_entry("DRBARMASS", 1000005);
   LOCALDRBAR(MSd)(3) = slha_io.read_entry("DRBARMASS", 2000001);
   LOCALDRBAR(MSd)(4) = slha_io.read_entry("DRBARMASS", 2000003);
   LOCALDRBAR(MSd)(5) = slha_io.read_entry("DRBARMASS", 2000005);
   LOCALDRBAR(MSv)(0) = slha_io.read_entry("DRBARMASS", 1000012);
   LOCALDRBAR(MSv)(1) = slha_io.read_entry("DRBARMASS", 1000014);
   LOCALDRBAR(MSv)(2) = slha_io.read_entry("DRBARMASS", 1000016);
   LOCALDRBAR(MSu)(0) = slha_io.read_entry("DRBARMASS", 1000002);
   LOCALDRBAR(MSu)(1) = slha_io.read_entry("DRBARMASS", 1000004);
   LOCALDRBAR(MSu)(2) = slha_io.read_entry("DRBARMASS", 1000006);
   LOCALDRBAR(MSu)(3) = slha_io.read_entry("DRBARMASS", 2000002);
   LOCALDRBAR(MSu)(4) = slha_io.read_entry("DRBARMASS", 2000004);
   LOCALDRBAR(MSu)(5) = slha_io.read_entry("DRBARMASS", 2000006);
   LOCALDRBAR(MSe)(0) = slha_io.read_entry("DRBARMASS", 1000011);
   LOCALDRBAR(MSe)(1) = slha_io.read_entry("DRBARMASS", 1000013);
   LOCALDRBAR(MSe)(2) = slha_io.read_entry("DRBARMASS", 1000015);
   LOCALDRBAR(MSe)(3) = slha_io.read_entry("DRBARMASS", 2000011);
   LOCALDRBAR(MSe)(4) = slha_io.read_entry("DRBARMASS", 2000013);
   LOCALDRBAR(MSe)(5) = slha_io.read_entry("DRBARMASS", 2000015);
   LOCALDRBAR(MSDX)(0) = slha_io.read_entry("DRBARMASS", 1000051);
   LOCALDRBAR(MSDX)(1) = slha_io.read_entry("DRBARMASS", 2000051);
   LOCALDRBAR(MSDX)(2) = slha_io.read_entry("DRBARMASS", 1000052);
   LOCALDRBAR(MSDX)(3) = slha_io.read_entry("DRBARMASS", 2000052);
   LOCALDRBAR(MSDX)(4) = slha_io.read_entry("DRBARMASS", 1000053);
   LOCALDRBAR(MSDX)(5) = slha_io.read_entry("DRBARMASS", 2000053);
   LOCALDRBAR(Mhh)(0) = slha_io.read_entry("DRBARMASS", 25);
   LOCALDRBAR(Mhh)(1) = slha_io.read_entry("DRBARMASS", 35);
   LOCALDRBAR(Mhh)(2) = slha_io.read_entry("DRBARMASS", 45);
   LOCALDRBAR(Mhh)(3) = slha_io.read_entry("DRBARMASS", 55);
   LOCALDRBAR(Mhh)(4) = slha_io.read_entry("DRBARMASS", 65);
   LOCALDRBAR(MAh)(2) = slha_io.read_entry("DRBARMASS", 91191138);
   LOCALDRBAR(MAh)(3) = slha_io.read_entry("DRBARMASS", 36);
   LOCALDRBAR(MAh)(4) = slha_io.read_entry("DRBARMASS", 91191137);
   LOCALDRBAR(MHpm)(1) = slha_io.read_entry("DRBARMASS", 37);
   LOCALDRBAR(MChi)(0) = slha_io.read_entry("DRBARMASS", 1000022);
   LOCALDRBAR(MChi)(1) = slha_io.read_entry("DRBARMASS", 1000023);
   LOCALDRBAR(MChi)(2) = slha_io.read_entry("DRBARMASS", 1000025);
   LOCALDRBAR(MChi)(3) = slha_io.read_entry("DRBARMASS", 1000035);
   LOCALDRBAR(MChi)(4) = slha_io.read_entry("DRBARMASS", 1000045);
   LOCALDRBAR(MChi)(5) = slha_io.read_entry("DRBARMASS", 1000055);
   LOCALDRBAR(MChi)(6) = slha_io.read_entry("DRBARMASS", 1000065);
   LOCALDRBAR(MChi)(7) = slha_io.read_entry("DRBARMASS", 1000075);
   LOCALDRBAR(MCha)(0) = slha_io.read_entry("DRBARMASS", 1000024);
   LOCALDRBAR(MCha)(1) = slha_io.read_entry("DRBARMASS", 1000037);
   LOCALDRBAR(MFe)(0) = slha_io.read_entry("DRBARMASS", 11);
   LOCALDRBAR(MFe)(1) = slha_io.read_entry("DRBARMASS", 13);
   LOCALDRBAR(MFe)(2) = slha_io.read_entry("DRBARMASS", 15);
   LOCALDRBAR(MFd)(0) = slha_io.read_entry("DRBARMASS", 1);
   LOCALDRBAR(MFd)(1) = slha_io.read_entry("DRBARMASS", 3);
   LOCALDRBAR(MFd)(2) = slha_io.read_entry("DRBARMASS", 5);
   LOCALDRBAR(MFu)(0) = slha_io.read_entry("DRBARMASS", 2);
   LOCALDRBAR(MFu)(1) = slha_io.read_entry("DRBARMASS", 4);
   LOCALDRBAR(MFu)(2) = slha_io.read_entry("DRBARMASS", 6);
   LOCALDRBAR(MFDX)(0) = slha_io.read_entry("DRBARMASS", 51);
   LOCALDRBAR(MFDX)(1) = slha_io.read_entry("DRBARMASS", 52);
   LOCALDRBAR(MFDX)(2) = slha_io.read_entry("DRBARMASS", 53);
   LOCALDRBAR(MSHI0)(0) = slha_io.read_entry("DRBARMASS", 82);
   LOCALDRBAR(MSHI0)(1) = slha_io.read_entry("DRBARMASS", 86);
   LOCALDRBAR(MSHI0)(2) = slha_io.read_entry("DRBARMASS", 84);
   LOCALDRBAR(MSHI0)(3) = slha_io.read_entry("DRBARMASS", 88);
   LOCALDRBAR(MSHI0)(4) = slha_io.read_entry("DRBARMASS", 9994453);
   LOCALDRBAR(MSHI0)(5) = slha_io.read_entry("DRBARMASS", 9994454);
   LOCALDRBAR(MSHI0)(6) = slha_io.read_entry("DRBARMASS", 9994455);
   LOCALDRBAR(MSHIPM)(0) = slha_io.read_entry("DRBARMASS", 81);
   LOCALDRBAR(MSHIPM)(1) = slha_io.read_entry("DRBARMASS", 85);
   LOCALDRBAR(MSHIPM)(2) = slha_io.read_entry("DRBARMASS", 83);
   LOCALDRBAR(MSHIPM)(3) = slha_io.read_entry("DRBARMASS", 87);
   LOCALDRBAR(MChaI)(0) = slha_io.read_entry("DRBARMASS", 1000088);
   LOCALDRBAR(MChaI)(1) = slha_io.read_entry("DRBARMASS", 1000089);
   LOCALDRBAR(MChiI)(0) = slha_io.read_entry("DRBARMASS", 1000081);
   LOCALDRBAR(MChiI)(1) = slha_io.read_entry("DRBARMASS", 1000082);
   LOCALDRBAR(MChiI)(2) = slha_io.read_entry("DRBARMASS", 1000083);
   LOCALDRBAR(MChiI)(3) = slha_io.read_entry("DRBARMASS", 1000084);
   LOCALDRBAR(MChiI)(4) = slha_io.read_entry("DRBARMASS", 1000085);
   LOCALDRBAR(MChiI)(5) = slha_io.read_entry("DRBARMASS", 1000086);
   LOCALDRBAR(MChiI)(6) = slha_io.read_entry("DRBARMASS", 1000087);
   LOCALDRBAR(MSHp0)(0) = slha_io.read_entry("DRBARMASS", 92);
   LOCALDRBAR(MSHp0)(1) = slha_io.read_entry("DRBARMASS", 94);
   LOCALDRBAR(MSHpp)(0) = slha_io.read_entry("DRBARMASS", 91);
   LOCALDRBAR(MSHpp)(1) = slha_io.read_entry("DRBARMASS", 93);
   LOCALDRBAR(MChiP)(0) = slha_io.read_entry("DRBARMASS", 1000092);
   LOCALDRBAR(MChiP)(1) = slha_io.read_entry("DRBARMASS", 1000094);
   LOCALDRBAR(MVWm) = slha_io.read_entry("DRBARMASS", 24);

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
void CSE6SSM_slha_io::fill_physical(CSE6SSM_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(ZD);
      slha_io.read_block("DSQMIX", ZD);
      LOCALPHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZU);
      slha_io.read_block("USQMIX", ZU);
      LOCALPHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZE);
      slha_io.read_block("SELMIX", ZE);
      LOCALPHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZV);
      slha_io.read_block("SNUMIX", ZV);
      LOCALPHYSICAL(ZV) = ZV;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDX);
      slha_io.read_block("ESIXZDX", ZDX);
      LOCALPHYSICAL(ZDX) = ZDX;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("NMHMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("NMAMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN);
      slha_io.read_block("NMNMIX", ZN);
      LOCALPHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZNp);
      slha_io.read_block("ZNPMIX", ZNp);
      LOCALPHYSICAL(ZNp) = ZNp;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZNI);
      slha_io.read_block("ZNIMIX", ZNI);
      LOCALPHYSICAL(ZNI) = ZNI;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZMI);
      slha_io.read_block("ESIXZMI", ZMI);
      LOCALPHYSICAL(ZMI) = ZMI;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZPI);
      slha_io.read_block("ESIXZPI", ZPI);
      LOCALPHYSICAL(ZPI) = ZPI;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      LOCALPHYSICAL(UP) = UP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      LOCALPHYSICAL(UM) = UM;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZEL);
      slha_io.read_block("UELMIX", ZEL);
      LOCALPHYSICAL(ZEL) = ZEL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZER);
      slha_io.read_block("UERMIX", ZER);
      LOCALPHYSICAL(ZER) = ZER;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDL);
      slha_io.read_block("UDLMIX", ZDL);
      LOCALPHYSICAL(ZDL) = ZDL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDR);
      slha_io.read_block("UDRMIX", ZDR);
      LOCALPHYSICAL(ZDR) = ZDR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZUL);
      slha_io.read_block("UULMIX", ZUL);
      LOCALPHYSICAL(ZUL) = ZUL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZUR);
      slha_io.read_block("UURMIX", ZUR);
      LOCALPHYSICAL(ZUR) = ZUR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDXL);
      slha_io.read_block("ESIXZXL", ZDXL);
      LOCALPHYSICAL(ZDXL) = ZDXL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDXR);
      slha_io.read_block("ESIXZXR", ZDXR);
      LOCALPHYSICAL(ZDXR) = ZDXR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UHp0);
      slha_io.read_block("UHNPMIX", UHp0);
      LOCALPHYSICAL(UHp0) = UHp0;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UHpp);
      slha_io.read_block("UHPPMIX", UHpp);
      LOCALPHYSICAL(UHpp) = UHpp;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UHI0);
      slha_io.read_block("UHNIMIX", UHI0);
      LOCALPHYSICAL(UHI0) = UHI0;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UHIPM);
      slha_io.read_block("UHPPMIX", UHIPM);
      LOCALPHYSICAL(UHIPM) = UHIPM;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MChaP) = slha_io.read_entry("MASS", 1000091);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(MVZp) = slha_io.read_entry("MASS", 31);
   LOCALPHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   LOCALPHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 1000003);
   LOCALPHYSICAL(MSd)(2) = slha_io.read_entry("MASS", 1000005);
   LOCALPHYSICAL(MSd)(3) = slha_io.read_entry("MASS", 2000001);
   LOCALPHYSICAL(MSd)(4) = slha_io.read_entry("MASS", 2000003);
   LOCALPHYSICAL(MSd)(5) = slha_io.read_entry("MASS", 2000005);
   LOCALPHYSICAL(MSv)(0) = slha_io.read_entry("MASS", 1000012);
   LOCALPHYSICAL(MSv)(1) = slha_io.read_entry("MASS", 1000014);
   LOCALPHYSICAL(MSv)(2) = slha_io.read_entry("MASS", 1000016);
   LOCALPHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   LOCALPHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 1000004);
   LOCALPHYSICAL(MSu)(2) = slha_io.read_entry("MASS", 1000006);
   LOCALPHYSICAL(MSu)(3) = slha_io.read_entry("MASS", 2000002);
   LOCALPHYSICAL(MSu)(4) = slha_io.read_entry("MASS", 2000004);
   LOCALPHYSICAL(MSu)(5) = slha_io.read_entry("MASS", 2000006);
   LOCALPHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   LOCALPHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 1000013);
   LOCALPHYSICAL(MSe)(2) = slha_io.read_entry("MASS", 1000015);
   LOCALPHYSICAL(MSe)(3) = slha_io.read_entry("MASS", 2000011);
   LOCALPHYSICAL(MSe)(4) = slha_io.read_entry("MASS", 2000013);
   LOCALPHYSICAL(MSe)(5) = slha_io.read_entry("MASS", 2000015);
   LOCALPHYSICAL(MSDX)(0) = slha_io.read_entry("MASS", 1000051);
   LOCALPHYSICAL(MSDX)(1) = slha_io.read_entry("MASS", 2000051);
   LOCALPHYSICAL(MSDX)(2) = slha_io.read_entry("MASS", 1000052);
   LOCALPHYSICAL(MSDX)(3) = slha_io.read_entry("MASS", 2000052);
   LOCALPHYSICAL(MSDX)(4) = slha_io.read_entry("MASS", 1000053);
   LOCALPHYSICAL(MSDX)(5) = slha_io.read_entry("MASS", 2000053);
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   LOCALPHYSICAL(Mhh)(3) = slha_io.read_entry("MASS", 55);
   LOCALPHYSICAL(Mhh)(4) = slha_io.read_entry("MASS", 65);
   LOCALPHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 91191138);
   LOCALPHYSICAL(MAh)(3) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MAh)(4) = slha_io.read_entry("MASS", 91191137);
   LOCALPHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   LOCALPHYSICAL(MChi)(4) = slha_io.read_entry("MASS", 1000045);
   LOCALPHYSICAL(MChi)(5) = slha_io.read_entry("MASS", 1000055);
   LOCALPHYSICAL(MChi)(6) = slha_io.read_entry("MASS", 1000065);
   LOCALPHYSICAL(MChi)(7) = slha_io.read_entry("MASS", 1000075);
   LOCALPHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   LOCALPHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MFDX)(0) = slha_io.read_entry("MASS", 51);
   LOCALPHYSICAL(MFDX)(1) = slha_io.read_entry("MASS", 52);
   LOCALPHYSICAL(MFDX)(2) = slha_io.read_entry("MASS", 53);
   LOCALPHYSICAL(MSHI0)(0) = slha_io.read_entry("MASS", 82);
   LOCALPHYSICAL(MSHI0)(1) = slha_io.read_entry("MASS", 86);
   LOCALPHYSICAL(MSHI0)(2) = slha_io.read_entry("MASS", 84);
   LOCALPHYSICAL(MSHI0)(3) = slha_io.read_entry("MASS", 88);
   LOCALPHYSICAL(MSHI0)(4) = slha_io.read_entry("MASS", 9994453);
   LOCALPHYSICAL(MSHI0)(5) = slha_io.read_entry("MASS", 9994454);
   LOCALPHYSICAL(MSHI0)(6) = slha_io.read_entry("MASS", 9994455);
   LOCALPHYSICAL(MSHIPM)(0) = slha_io.read_entry("MASS", 81);
   LOCALPHYSICAL(MSHIPM)(1) = slha_io.read_entry("MASS", 85);
   LOCALPHYSICAL(MSHIPM)(2) = slha_io.read_entry("MASS", 83);
   LOCALPHYSICAL(MSHIPM)(3) = slha_io.read_entry("MASS", 87);
   LOCALPHYSICAL(MChaI)(0) = slha_io.read_entry("MASS", 1000088);
   LOCALPHYSICAL(MChaI)(1) = slha_io.read_entry("MASS", 1000089);
   LOCALPHYSICAL(MChiI)(0) = slha_io.read_entry("MASS", 1000081);
   LOCALPHYSICAL(MChiI)(1) = slha_io.read_entry("MASS", 1000082);
   LOCALPHYSICAL(MChiI)(2) = slha_io.read_entry("MASS", 1000083);
   LOCALPHYSICAL(MChiI)(3) = slha_io.read_entry("MASS", 1000084);
   LOCALPHYSICAL(MChiI)(4) = slha_io.read_entry("MASS", 1000085);
   LOCALPHYSICAL(MChiI)(5) = slha_io.read_entry("MASS", 1000086);
   LOCALPHYSICAL(MChiI)(6) = slha_io.read_entry("MASS", 1000087);
   LOCALPHYSICAL(MSHp0)(0) = slha_io.read_entry("MASS", 92);
   LOCALPHYSICAL(MSHp0)(1) = slha_io.read_entry("MASS", 94);
   LOCALPHYSICAL(MSHpp)(0) = slha_io.read_entry("MASS", 91);
   LOCALPHYSICAL(MSHpp)(1) = slha_io.read_entry("MASS", 93);
   LOCALPHYSICAL(MChiP)(0) = slha_io.read_entry("MASS", 1000092);
   LOCALPHYSICAL(MChiP)(1) = slha_io.read_entry("MASS", 1000094);
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double CSE6SSM_slha_io::read_scale() const
{
   double scale = 0.;

   for (unsigned i = 0; i < NUMBER_OF_DRBAR_BLOCKS; i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
