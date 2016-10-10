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

// File generated at Fri 10 Jul 2015 12:04:45

#include "CMSSM_slha_io.hpp"
#include "CMSSM_two_scale_input_parameters.hpp"
#include "CMSSM_semi_two_scale_input_parameters.hpp"
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

char const * const CMSSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "Te", "Td", "Tu", "HMIX", "MSQ2", "MSE2",
   "MSL2", "MSU2", "MSD2", "MSOFT" }
;

CMSSM_slha_io::CMSSM_slha_io()
   : slha_io()
{
}

void CMSSM_slha_io::clear()
{
   slha_io.clear();
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CMSSM_slha_io::set_extpar(const CMSSM_input_parameters<Two_scale>& input)
{

}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CMSSM_slha_io::set_extpar(const CMSSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   if (input.MuInput_at_MS) {
      extpar << FORMAT_ELEMENT(0, -1, "Mu input at SUSY scale");
   }
   extpar << FORMAT_ELEMENT(23, input.MuInput, "MuInput");
   slha_io.set_block(extpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CMSSM_slha_io::set_minpar(const CMSSM_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.m0, "m0");
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(4, input.SignMu, "SignMu");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CMSSM_slha_io::set_minpar(const CMSSM_semianalytic_input_parameters<Two_scale>& input)
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
void CMSSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void CMSSM_slha_io::set_spinfo(const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems)
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
 * Stores the model DR-bar masses in the SLHA object.
 *
 * @param drbar struct of DR-bar masses and mixings
 *
 * @param scale scale at which the DR-bar masses are evaluated
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void CMSSM_slha_io::set_drbar_mass(const MSSM_physical& drbar,
                                   double scale,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block DRBARMASS Q= " << FORMAT_SCALE(scale) << "\n"
      << FORMAT_MASS(1000021, LOCALDRBAR(MGlu), "DRbarGlu")
      << FORMAT_MASS(24, LOCALDRBAR(MVWm), "DRbarVWm")
      << FORMAT_MASS(1000024, LOCALDRBAR(MCha(0)), "DRbarCha(1)")
      << FORMAT_MASS(1000037, LOCALDRBAR(MCha(1)), "DRbarCha(2)")
      << FORMAT_MASS(25, LOCALDRBAR(Mhh(0)), "DRbarhh(1)")
      << FORMAT_MASS(35, LOCALDRBAR(Mhh(1)), "DRbarhh(2)")
      << FORMAT_MASS(37, LOCALDRBAR(MHpm(1)), "DRbarHpm(2)")
      << FORMAT_MASS(36, LOCALDRBAR(MAh(1)), "DRbarAh(2)")
      << FORMAT_MASS(1000012, LOCALDRBAR(MSv(0)), "DRbarSv(1)")
      << FORMAT_MASS(1000014, LOCALDRBAR(MSv(1)), "DRbarSv(2)")
      << FORMAT_MASS(1000016, LOCALDRBAR(MSv(2)), "DRbarSv(3)")
      << FORMAT_MASS(1000022, LOCALDRBAR(MChi(0)), "DRbarChi(1)")
      << FORMAT_MASS(1000023, LOCALDRBAR(MChi(1)), "DRbarChi(2)")
      << FORMAT_MASS(1000025, LOCALDRBAR(MChi(2)), "DRbarChi(3)")
      << FORMAT_MASS(1000035, LOCALDRBAR(MChi(3)), "DRbarChi(4)")
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
 * @param drbar struct of DR-bar parameters
 *
 * @param scale scale at which the DR-bar mixings are evaluated
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void CMSSM_slha_io::set_drbar_mixing_matrices(const MSSM_physical& drbar,
                                              double scale,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("DRBARUMIX", LOCALDRBAR(UM), "DRbarUM", scale);
   slha_io.set_block("DRBARVMIX", LOCALDRBAR(UP), "DRbarUP", scale);
   slha_io.set_block("DRBARPSEUDOSCALARMIX", LOCALDRBAR(ZA), "DRbarZA", scale);
   slha_io.set_block("DRBARDSQMIX", LOCALDRBAR(ZD), "DRbarZD", scale);
   slha_io.set_block("DRBARSELMIX", LOCALDRBAR(ZE), "DRbarZE", scale);
   slha_io.set_block("DRBARSCALARMIX", LOCALDRBAR(ZH), "DRbarZH", scale);
   slha_io.set_block("DRBARNMIX", LOCALDRBAR(ZN), "DRbarZN", scale);
   slha_io.set_block("DRBARCHARGEMIX", LOCALDRBAR(ZP), "DRbarZP", scale);
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
void CMSSM_slha_io::set_mass(const MSSM_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha(2)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHpm(1)), "Hpm(2)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(1)), "Ah(2)")
      << FORMAT_MASS(1000012, LOCALPHYSICAL(MSv(0)), "Sv(1)")
      << FORMAT_MASS(1000014, LOCALPHYSICAL(MSv(1)), "Sv(2)")
      << FORMAT_MASS(1000016, LOCALPHYSICAL(MSv(2)), "Sv(3)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
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
void CMSSM_slha_io::set_mixing_matrices(const MSSM_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("PSEUDOSCALARMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", LOCALPHYSICAL(ZD), "ZD");
   slha_io.set_block("SELMIX", LOCALPHYSICAL(ZE), "ZE");
   slha_io.set_block("SCALARMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("NMIX", LOCALPHYSICAL(ZN), "ZN");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");
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

void CMSSM_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void CMSSM_slha_io::set_pmns(
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
void CMSSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double CMSSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void CMSSM_slha_io::read_from_file(const std::string& file_name)
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
void CMSSM_slha_io::fill(CMSSM_input_parameters<Two_scale>& input) const
{
   void (*fill_two_scale_minpar_tuple) (CMSSM_input_parameters<Two_scale>&,
                                        int, double) = &CMSSM_slha_io::fill_minpar_tuple;
   void (*fill_two_scale_extpar_tuple) (CMSSM_input_parameters<Two_scale>&,
                                        int, double) = &CMSSM_slha_io::fill_extpar_tuple;

   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_two_scale_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_two_scale_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);


}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void CMSSM_slha_io::fill(CMSSM_semianalytic_input_parameters<Two_scale>& input) const
{
   void (*fill_semi_two_scale_minpar_tuple) (CMSSM_semianalytic_input_parameters<Two_scale>&,
                                             int, double) = &CMSSM_slha_io::fill_minpar_tuple;
   void (*fill_semi_two_scale_extpar_tuple) (CMSSM_semianalytic_input_parameters<Two_scale>&,
                                             int, double) = &CMSSM_slha_io::fill_extpar_tuple;

   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_semi_two_scale_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_semi_two_scale_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);


}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
void CMSSM_slha_io::fill_drbar_parameters(MSSM_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
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
   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_BMu(slha_io.read_entry("HMIX", 101));
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
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 */
void CMSSM_slha_io::fill(MSSM_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   MSSM_physical drbar_hk;
   fill_drbar(drbar_hk);
   drbar_hk.convert_to_hk();
   model.get_drbar_masses() = drbar_hk;

   MSSM_physical physical_hk;
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
void CMSSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&CMSSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void CMSSM_slha_io::fill_minpar_tuple(CMSSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.m0 = value; break;
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 4: input.SignMu = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CMSSM_slha_io::fill_minpar_tuple(CMSSM_semianalytic_input_parameters<Two_scale>& input,
                                      int key, double value)
{
   switch (key) {
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CMSSM_slha_io::fill_extpar_tuple(CMSSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CMSSM_slha_io::fill_extpar_tuple(CMSSM_semianalytic_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 0: {
      if (value < 0) {
         input.MuInput_at_MS = true;
      } else {
         // N.B. this is not compliant with SLHA2
         input.MuInput_at_MS = false;
      }
      break;
   }
   case 23: input.MuInput = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CMSSM_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
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
void CMSSM_slha_io::fill_drbar(MSSM_physical& drbar) const
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
      DEFINE_DRBAR_PARAMETER(ZH);
      slha_io.read_block("DRBARSCALARMIX", ZH);
      LOCALDRBAR(ZH) = ZH;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZA);
      slha_io.read_block("DRBARPSEUDOSCALARMIX", ZA);
      LOCALDRBAR(ZA) = ZA;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZP);
      slha_io.read_block("DRBARCHARGEMIX", ZP);
      LOCALDRBAR(ZP) = ZP;
   }
   {
      DEFINE_DRBAR_PARAMETER(ZN);
      slha_io.read_block("DRBARNMIX", ZN);
      LOCALDRBAR(ZN) = ZN;
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

   LOCALDRBAR(MVG) = slha_io.read_entry("DRBARMASS", 21);
   LOCALDRBAR(MGlu) = slha_io.read_entry("DRBARMASS", 1000021);
   LOCALDRBAR(MFv)(0) = slha_io.read_entry("DRBARMASS", 12);
   LOCALDRBAR(MFv)(1) = slha_io.read_entry("DRBARMASS", 14);
   LOCALDRBAR(MFv)(2) = slha_io.read_entry("DRBARMASS", 16);
   LOCALDRBAR(MVP) = slha_io.read_entry("DRBARMASS", 22);
   LOCALDRBAR(MVZ) = slha_io.read_entry("DRBARMASS", 23);
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
   LOCALDRBAR(Mhh)(0) = slha_io.read_entry("DRBARMASS", 25);
   LOCALDRBAR(Mhh)(1) = slha_io.read_entry("DRBARMASS", 35);
   LOCALDRBAR(MAh)(1) = slha_io.read_entry("DRBARMASS", 36);
   LOCALDRBAR(MHpm)(1) = slha_io.read_entry("DRBARMASS", 37);
   LOCALDRBAR(MChi)(0) = slha_io.read_entry("DRBARMASS", 1000022);
   LOCALDRBAR(MChi)(1) = slha_io.read_entry("DRBARMASS", 1000023);
   LOCALDRBAR(MChi)(2) = slha_io.read_entry("DRBARMASS", 1000025);
   LOCALDRBAR(MChi)(3) = slha_io.read_entry("DRBARMASS", 1000035);
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
   LOCALDRBAR(MVWm) = slha_io.read_entry("DRBARMASS", 24);

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
void CMSSM_slha_io::fill_physical(MSSM_physical& physical) const
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
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN);
      slha_io.read_block("NMIX", ZN);
      LOCALPHYSICAL(ZN) = ZN;
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

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
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
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
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
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double CMSSM_slha_io::read_scale() const
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
