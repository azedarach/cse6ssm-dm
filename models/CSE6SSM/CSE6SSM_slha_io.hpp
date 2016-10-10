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

#ifndef CSE6SSM_SLHA_IO_H
#define CSE6SSM_SLHA_IO_H

#include "CSE6SSM_semi_two_scale_model_slha.hpp"
#include "CSE6SSM_two_scale_model_slha.hpp"
#include "CSE6SSM_info.hpp"
#include "CSE6SSM_physical.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define SEMIANALYTICCOEFF(p,q) model.get_##p##_coeff_##q()
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

template <class T>
struct CSE6SSM_input_parameters;

template <class T>
struct CSE6SSM_semianalytic_input_parameters;

class Two_scale;
class Spectrum_generator_settings;

struct CSE6SSM_scales {
   CSE6SSM_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class CSE6SSM_slha_io {
public:
   CSE6SSM_slha_io();
   ~CSE6SSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(CSE6SSM_input_parameters<Two_scale>&) const;
   void fill(CSE6SSM_semianalytic_input_parameters<Two_scale>&) const;
   void fill(CSE6SSM_mass_eigenstates&) const;
   template <class T> void fill(CSE6SSM_slha<T>&) const;
   template <class T> void fill(CSE6SSM_semianalytic_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const CSE6SSM_input_parameters<Two_scale>&);
   void set_extpar(const CSE6SSM_semianalytic_input_parameters<Two_scale>&);
   template <class T> void set_extra(const CSE6SSM_slha<T>&, const CSE6SSM_scales&);
   template <class T> void set_extra(const CSE6SSM_semianalytic_slha<T>&, const CSE6SSM_scales&);
   template <class T> void set_extra(const CSE6SSM_semianalytic_slha<T>&, const CSE6SSM_scales&, double MhhEFT);
   void set_minpar(const CSE6SSM_input_parameters<Two_scale>&);
   void set_minpar(const CSE6SSM_semianalytic_input_parameters<Two_scale>&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const CSE6SSM_slha<T>&);
   template <class T> void set_spectrum(const CSE6SSM<T>&);
   template <class T> void set_spectrum(const CSE6SSM_semianalytic_slha<T>&);
   template <class T> void set_spectrum(const CSE6SSM_semianalytic<T>&);
   void set_spinfo(const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(CSE6SSM_input_parameters<Two_scale>&, int, double);
   static void fill_extpar_tuple(CSE6SSM_input_parameters<Two_scale>&, int, double);
   static void fill_minpar_tuple(CSE6SSM_semianalytic_input_parameters<Two_scale>&, int, double);
   static void fill_extpar_tuple(CSE6SSM_semianalytic_input_parameters<Two_scale>&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CSE6SSM_slha<T>&, const QedQcd&, const CSE6SSM_scales&);
   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CSE6SSM_semianalytic_slha<T>&, const QedQcd&, const CSE6SSM_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CSE6SSM_slha<T>&, const QedQcd&);
   template <class T>
   static SLHAea::Coll fill_slhaea(const CSE6SSM_semianalytic_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CSE6SSM_slha<T>&, const QedQcd&, const CSE6SSM_scales&);
   template <class T>
   static SLHAea::Coll fill_slhaea(const CSE6SSM_semianalytic_slha<T>&, const QedQcd&, const CSE6SSM_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 59;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_drbar_mass(const CSE6SSM_physical&, double, bool);
   void set_drbar_mixing_matrices(const CSE6SSM_physical&, double, bool);
   void set_mass(const CSE6SSM_physical&, bool);
   void set_mixing_matrices(const CSE6SSM_physical&, bool);
   template <class T> void set_model_parameters(const CSE6SSM_slha<T>&);
   template <class T> void set_model_parameters(const CSE6SSM_semianalytic_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(CSE6SSM_mass_eigenstates&) const;
   void fill_drbar(CSE6SSM_physical&) const;
   void fill_physical(CSE6SSM_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CSE6SSM_slha_io::fill(CSE6SSM_slha<T>& model) const
{
   fill(static_cast<CSE6SSM_mass_eigenstates&>(model));
   fill_drbar(model.get_drbar_slha());
   fill_physical(model.get_physical_slha());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CSE6SSM_slha_io::fill(CSE6SSM_semianalytic_slha<T>& model) const
{
   fill(static_cast<CSE6SSM_mass_eigenstates&>(model));
   fill_drbar(model.get_drbar_slha());
   fill_physical(model.get_physical_slha());
}

template <class T>
void CSE6SSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CSE6SSM_slha<T>& model,
   const QedQcd& qedqcd, const CSE6SSM_scales& scales)
{
   CSE6SSM_slha_io slha_io;
   const CSE6SSM_input_parameters<T>& input = model.get_input();
   const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
void CSE6SSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CSE6SSM_semianalytic_slha<T>& model,
   const QedQcd& qedqcd, const CSE6SSM_scales& scales)
{
   CSE6SSM_slha_io slha_io;
   const CSE6SSM_semianalytic_input_parameters<T>& input = model.get_input();
   const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
SLHAea::Coll CSE6SSM_slha_io::fill_slhaea(
   const CSE6SSM_slha<T>& model, const QedQcd& qedqcd)
{
   CSE6SSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CSE6SSM_slha_io::fill_slhaea(
   const CSE6SSM_semianalytic_slha<T>& model, const QedQcd& qedqcd)
{
   CSE6SSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CSE6SSM_slha_io::fill_slhaea(
   const CSE6SSM_slha<T>& model, const QedQcd& qedqcd,
   const CSE6SSM_scales& scales)
{
   SLHAea::Coll slhaea;
   CSE6SSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

template <class T>
SLHAea::Coll CSE6SSM_slha_io::fill_slhaea(
   const CSE6SSM_semianalytic_slha<T>& model, const QedQcd& qedqcd,
   const CSE6SSM_scales& scales)
{
   SLHAea::Coll slhaea;
   CSE6SSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CSE6SSM_slha_io::set_model_parameters(const CSE6SSM_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(g1p)), "g1p")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block NCharge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(QS)), "QS")
         ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block("Te", MODELPARAMETER(TYe_slha), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd_slha), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu_slha), "TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(10, ArcTan((MODELPARAMETER(vu)) / (MODELPARAMETER(vd))), "Beta")
            << FORMAT_ELEMENT(30, (MODELPARAMETER(XiF)), "XiF")
            << FORMAT_ELEMENT(31, (MODELPARAMETER(MuPhi)), "MuPhi")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(KappaPr)), "KappaPr")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Sigmax)), "Sigmax")
            << FORMAT_ELEMENT(34, (MODELPARAMETER(LXiF)), "LXiF")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXHEYUK", MODELPARAMETER(hE), "hE", model.get_scale());
   {
      std::ostringstream block;
      block << "Block ESIXRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(42, (MODELPARAMETER(SigmaL)), "SigmaL")
            << FORMAT_ELEMENT(28, (MODELPARAMETER(TKappaPr)), "TKappaPr")
            << FORMAT_ELEMENT(29, (MODELPARAMETER(TSigmax)), "TSigmax")
            << FORMAT_ELEMENT(43, (MODELPARAMETER(TSigmaL)), "TSigmaL")
            << FORMAT_ELEMENT(30, (MODELPARAMETER(BMuPhi)), "BMuPhi")
            << FORMAT_ELEMENT(11, (MODELPARAMETER(vs)), "vs")
            << FORMAT_ELEMENT(12, (MODELPARAMETER(vsb)), "vsb")
            << FORMAT_ELEMENT(13, (MODELPARAMETER(vphi)), "vphi")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(0, (MODELPARAMETER(MuPr)), "MuPr")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMuPr)), "BMuPr")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXGDYUK", MODELPARAMETER(gD), "gD", model.get_scale());
   slha_io.set_block("ESIXFUYUK", MODELPARAMETER(fu), "fu", model.get_scale());
   slha_io.set_block("ESIXFDYUK", MODELPARAMETER(fd), "fd", model.get_scale());
   slha_io.set_block("ESIXTHETRI", MODELPARAMETER(ThE), "ThE", model.get_scale());
   slha_io.set_block("ESIXTGDTRI", MODELPARAMETER(TgD), "TgD", model.get_scale());
   slha_io.set_block("ESIXTFUTRI", MODELPARAMETER(Tfu), "Tfu", model.get_scale());
   slha_io.set_block("ESIXTFDTRI", MODELPARAMETER(Tfd), "Tfd", model.get_scale());
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2_slha), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2_slha), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2_slha), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2_slha), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2_slha), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(23, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(24, (MODELPARAMETER(msbar2)), "msbar2")
            << FORMAT_ELEMENT(25, (MODELPARAMETER(mphi2)), "mphi2")
            << FORMAT_ELEMENT(26, (MODELPARAMETER(mHp2)), "mHp2")
            << FORMAT_ELEMENT(27, (MODELPARAMETER(mHpbar2)), "mHpbar2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(MassBp)), "MassBp")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mX2", MODELPARAMETER(mDx2), "mDx2", model.get_scale());
   slha_io.set_block("mXBar2", MODELPARAMETER(mDxbar2), "mDxbar2", model.get_scale());
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXTKAPPA", MODELPARAMETER(TKappa), "TKappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());
   slha_io.set_block("ESIXTLAMBDA", MODELPARAMETER(TLambda12), "TLambda12", model.get_scale());
   {
      std::ostringstream block;
      block << "Block PHASES Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, Re((MODELPARAMETER(PhaseGlu))), "Re(PhaseGlu)")
         ;
      slha_io.set_block(block);
   }
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CSE6SSM_slha_io::set_model_parameters(const CSE6SSM_semianalytic_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(g1p)), "g1p")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block NCharge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(QS)), "QS")
         ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block("Te", MODELPARAMETER(TYe_slha), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd_slha), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu_slha), "TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(10, ArcTan((MODELPARAMETER(vu)) / (MODELPARAMETER(vd))), "Beta")
            << FORMAT_ELEMENT(30, (MODELPARAMETER(XiF)), "XiF")
            << FORMAT_ELEMENT(31, (MODELPARAMETER(MuPhi)), "MuPhi")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(KappaPr)), "KappaPr")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Sigmax)), "Sigmax")
            << FORMAT_ELEMENT(34, (MODELPARAMETER(LXiF)), "LXiF")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXHEYUK", MODELPARAMETER(hE), "hE", model.get_scale());
   {
      std::ostringstream block;
      block << "Block ESIXRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(42, (MODELPARAMETER(SigmaL)), "SigmaL")
            << FORMAT_ELEMENT(28, (MODELPARAMETER(TKappaPr)), "TKappaPr")
            << FORMAT_ELEMENT(29, (MODELPARAMETER(TSigmax)), "TSigmax")
            << FORMAT_ELEMENT(43, (MODELPARAMETER(TSigmaL)), "TSigmaL")
            << FORMAT_ELEMENT(30, (MODELPARAMETER(BMuPhi)), "BMuPhi")
            << FORMAT_ELEMENT(11, (MODELPARAMETER(vs)), "vs")
            << FORMAT_ELEMENT(12, (MODELPARAMETER(vsb)), "vsb")
            << FORMAT_ELEMENT(13, (MODELPARAMETER(vphi)), "vphi")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(0, (MODELPARAMETER(MuPr)), "MuPr")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMuPr)), "BMuPr")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXGDYUK", MODELPARAMETER(gD), "gD", model.get_scale());
   slha_io.set_block("ESIXFUYUK", MODELPARAMETER(fu), "fu", model.get_scale());
   slha_io.set_block("ESIXFDYUK", MODELPARAMETER(fd), "fd", model.get_scale());
   slha_io.set_block("ESIXTHETRI", MODELPARAMETER(ThE), "ThE", model.get_scale());
   slha_io.set_block("ESIXTGDTRI", MODELPARAMETER(TgD), "TgD", model.get_scale());
   slha_io.set_block("ESIXTFUTRI", MODELPARAMETER(Tfu), "Tfu", model.get_scale());
   slha_io.set_block("ESIXTFDTRI", MODELPARAMETER(Tfd), "Tfd", model.get_scale());
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2_slha), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2_slha), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2_slha), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2_slha), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2_slha), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(23, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(24, (MODELPARAMETER(msbar2)), "msbar2")
            << FORMAT_ELEMENT(25, (MODELPARAMETER(mphi2)), "mphi2")
            << FORMAT_ELEMENT(26, (MODELPARAMETER(mHp2)), "mHp2")
            << FORMAT_ELEMENT(27, (MODELPARAMETER(mHpbar2)), "mHpbar2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(MassBp)), "MassBp")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mX2", MODELPARAMETER(mDx2), "mDx2", model.get_scale());
   slha_io.set_block("mXBar2", MODELPARAMETER(mDxbar2), "mDxbar2", model.get_scale());
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXTKAPPA", MODELPARAMETER(TKappa), "TKappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());
   slha_io.set_block("ESIXTLAMBDA", MODELPARAMETER(TLambda12), "TLambda12", model.get_scale());
   {
      std::ostringstream block;
      block << "Block PHASES Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, Re((MODELPARAMETER(PhaseGlu))), "Re(PhaseGlu)")
         ;
      slha_io.set_block(block);
   }
}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CSE6SSM_slha_io::set_extra(
   const CSE6SSM_slha<T>& model, const CSE6SSM_scales& scales)
{
   const CSE6SSM_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(1, (SCALES(SUSYScale)), "SUSYScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EWSBOutputParameters Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, model.get_ewsb_output_parameter(0), "TanTheta")
            << FORMAT_ELEMENT(1, model.get_ewsb_output_parameter(1), "Lambdax")
            << FORMAT_ELEMENT(2, model.get_ewsb_output_parameter(2), "vphi")
            << FORMAT_ELEMENT(3, model.get_ewsb_output_parameter(3), "XiF")
            << FORMAT_ELEMENT(4, model.get_ewsb_output_parameter(4), "LXiF")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(2,2)/Yu(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(3,3)/Yu(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(2,2)/Yd(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(3,3)/Yd(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(2,2)/Ye(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(3,3)/Ye(3,3)")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CSE6SSM_slha_io::set_extra(
   const CSE6SSM_semianalytic_slha<T>& model, const CSE6SSM_scales& scales)
{
   const CSE6SSM_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(1, (SCALES(SUSYScale)), "SUSYScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EWSBOutputParameters Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, model.get_ewsb_output_parameter(0), "m0Sqr")
            << FORMAT_ELEMENT(1, model.get_ewsb_output_parameter(1), "TanTheta")
            << FORMAT_ELEMENT(2, model.get_ewsb_output_parameter(2), "vphi")
            << FORMAT_ELEMENT(3, model.get_ewsb_output_parameter(3), "XiF")
            << FORMAT_ELEMENT(4, model.get_ewsb_output_parameter(4), "LXiF")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(2,2)/Yu(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(3,3)/Yu(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(2,2)/Yd(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(3,3)/Yd(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(2,2)/Ye(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(3,3)/Ye(3,3)")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("TeAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYe), "Azero_coeff_TYe", model.get_scale());
   slha_io.set_block("Tem12Coeffs", SEMIANALYTICCOEFF(m12,TYe), "m12_coeff_TYe", model.get_scale());
   slha_io.set_block("TdAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYd), "Azero_coeff_TYd", model.get_scale());
   slha_io.set_block("Tdm12Coeffs", SEMIANALYTICCOEFF(m12,TYd), "m12_coeff_TYd", model.get_scale());
   slha_io.set_block("TuAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYu), "Azero_coeff_TYu", model.get_scale());
   slha_io.set_block("Tum12Coeffs", SEMIANALYTICCOEFF(m12,TYu), "m12_coeff_TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block ESIXRUNCoeffs Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(281, (SEMIANALYTICCOEFF(Azero,TKappaPr)), "Azero_coeff_TKappaPr") << '\n'
            << FORMAT_ELEMENT(282, (SEMIANALYTICCOEFF(m12,TKappaPr)), "m12_coeff_TKappaPr") << '\n'
            << FORMAT_ELEMENT(291, (SEMIANALYTICCOEFF(Azero,TSigmax)), "Azero_coeff_TSigmax") << '\n'
            << FORMAT_ELEMENT(292, (SEMIANALYTICCOEFF(m12,TSigmax)), "m12_coeff_TSigmax") << '\n'
            << FORMAT_ELEMENT(431, (SEMIANALYTICCOEFF(Azero,TSigmaL)), "Azero_coeff_TSigmaL") << '\n'
            << FORMAT_ELEMENT(432, (SEMIANALYTICCOEFF(m12,TSigmaL)), "m12_coeff_TSigmaL") << '\n'
            << FORMAT_ELEMENT(301, (SEMIANALYTICCOEFF(Azero,BMuPhi)), "Azero_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(302, (SEMIANALYTICCOEFF(m12,BMuPhi)), "m12_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(303, (SEMIANALYTICCOEFF(BMuPhi,BMuPhi)), "BMuPhi_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(304, (SEMIANALYTICCOEFF(BMuPr,BMuPhi)), "BMuPr_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(21, (SEMIANALYTICCOEFF(Azero,TLambdax)), "Azero_coeff_TLambdax") << '\n'
            << FORMAT_ELEMENT(22, (SEMIANALYTICCOEFF(m12,TLambdax)), "m12_coeff_TLambdax") << '\n'
            << FORMAT_ELEMENT(1011, (SEMIANALYTICCOEFF(Azero,BMuPr)), "Azero_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1012, (SEMIANALYTICCOEFF(m12,BMuPr)), "m12_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1013, (SEMIANALYTICCOEFF(BMuPhi,BMuPr)), "BMuPhi_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1014, (SEMIANALYTICCOEFF(BMuPr,BMuPr)), "BMuPr_coeff_BMuPr") << '\n'
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ThEAzeroCoeffs", SEMIANALYTICCOEFF(Azero,ThE), "Azero_coeff_ThE", model.get_scale());
   slha_io.set_block("ThEm12Coeffs", SEMIANALYTICCOEFF(m12,ThE), "m12_coeff_ThE", model.get_scale());
   slha_io.set_block("TgDAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TgD), "Azero_coeff_TgD", model.get_scale());
   slha_io.set_block("TgDm12Coeffs", SEMIANALYTICCOEFF(m12,TgD), "m12_coeff_TgD", model.get_scale());
   slha_io.set_block("TfuAzeroCoeffs", SEMIANALYTICCOEFF(Azero,Tfu), "Azero_coeff_Tfu", model.get_scale());
   slha_io.set_block("Tfum12Coeffs", SEMIANALYTICCOEFF(m12,Tfu), "m12_coeff_Tfu", model.get_scale());
   slha_io.set_block("TfdAzeroCoeffs", SEMIANALYTICCOEFF(Azero,Tfd), "Azero_coeff_Tfd", model.get_scale());
   slha_io.set_block("Tfdm12Coeffs", SEMIANALYTICCOEFF(m12,Tfd), "m12_coeff_Tfd", model.get_scale());
   slha_io.set_block("mq2m02Coeffs", SEMIANALYTICCOEFF(m02,mq2), "m02_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2m122Coeffs", SEMIANALYTICCOEFF(m122,mq2), "m122_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mq2), "Azerom12_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mq2), "Azero2_coeff_mq2", model.get_scale());
   slha_io.set_block("me2m02Coeffs", SEMIANALYTICCOEFF(m02,me2), "m02_coeff_me2", model.get_scale());
   slha_io.set_block("me2m122Coeffs", SEMIANALYTICCOEFF(m122,me2), "m122_coeff_me2", model.get_scale());
   slha_io.set_block("me2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,me2), "Azerom12_coeff_me2", model.get_scale());
   slha_io.set_block("me2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,me2), "Azero2_coeff_me2", model.get_scale());
   slha_io.set_block("ml2m02Coeffs", SEMIANALYTICCOEFF(m02,ml2), "m02_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2m122Coeffs", SEMIANALYTICCOEFF(m122,ml2), "m122_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,ml2), "Azerom12_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,ml2), "Azero2_coeff_ml2", model.get_scale());
   slha_io.set_block("mu2m02Coeffs", SEMIANALYTICCOEFF(m02,mu2), "m02_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2m122Coeffs", SEMIANALYTICCOEFF(m122,mu2), "m122_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mu2), "Azerom12_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mu2), "Azero2_coeff_mu2", model.get_scale());
   slha_io.set_block("md2m02Coeffs", SEMIANALYTICCOEFF(m02,md2), "m02_coeff_md2", model.get_scale());
   slha_io.set_block("md2m122Coeffs", SEMIANALYTICCOEFF(m122,md2), "m122_coeff_md2", model.get_scale());
   slha_io.set_block("md2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,md2), "Azerom12_coeff_md2", model.get_scale());
   slha_io.set_block("md2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,md2), "Azero2_coeff_md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFTCoeffs Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(211, (SEMIANALYTICCOEFF(m02,mHd2)), "m02_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(212, (SEMIANALYTICCOEFF(m122,mHd2)), "m122_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(213, (SEMIANALYTICCOEFF(Azerom12,mHd2)), "Azerom12_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(214, (SEMIANALYTICCOEFF(Azero2,mHd2)), "Azero2_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(221, (SEMIANALYTICCOEFF(m02,mHu2)), "m02_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(222, (SEMIANALYTICCOEFF(m122,mHu2)), "m122_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(223, (SEMIANALYTICCOEFF(Azerom12,mHu2)), "Azerom12_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(224, (SEMIANALYTICCOEFF(Azero2,mHu2)), "Azero2_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(231, (SEMIANALYTICCOEFF(m02,ms2)), "m02_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(232, (SEMIANALYTICCOEFF(m122,ms2)), "m122_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(233, (SEMIANALYTICCOEFF(Azerom12,ms2)), "Azerom12_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(234, (SEMIANALYTICCOEFF(Azero2,ms2)), "Azero2_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(241, (SEMIANALYTICCOEFF(m02,msbar2)), "m02_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(242, (SEMIANALYTICCOEFF(m122,msbar2)), "m122_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(243, (SEMIANALYTICCOEFF(Azerom12,msbar2)), "Azerom12_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(244, (SEMIANALYTICCOEFF(Azero2,msbar2)), "Azero2_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(251, (SEMIANALYTICCOEFF(m02,mphi2)), "m02_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(252, (SEMIANALYTICCOEFF(m122,mphi2)), "m122_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(253, (SEMIANALYTICCOEFF(Azerom12,mphi2)), "Azerom12_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(254, (SEMIANALYTICCOEFF(Azero2,mphi2)), "Azero2_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(261, (SEMIANALYTICCOEFF(m02,mHp2)), "m02_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(262, (SEMIANALYTICCOEFF(m122,mHp2)), "m122_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(263, (SEMIANALYTICCOEFF(Azerom12,mHp2)), "Azerom12_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(264, (SEMIANALYTICCOEFF(Azero2,mHp2)), "Azero2_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(271, (SEMIANALYTICCOEFF(m02,mHpbar2)), "m02_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(272, (SEMIANALYTICCOEFF(m122,mHpbar2)), "m122_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(273, (SEMIANALYTICCOEFF(Azerom12,mHpbar2)), "Azerom12_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(274, (SEMIANALYTICCOEFF(Azero2,mHpbar2)), "Azero2_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(11, (SEMIANALYTICCOEFF(Azero,MassB)), "Azero_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(12, (SEMIANALYTICCOEFF(m12,MassB)), "m12_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(21, (SEMIANALYTICCOEFF(Azero,MassWB)), "Azero_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(22, (SEMIANALYTICCOEFF(m12,MassWB)), "m12_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(31, (SEMIANALYTICCOEFF(Azero,MassG)), "Azero_coeff_MassG") << '\n'
            << FORMAT_ELEMENT(32, (SEMIANALYTICCOEFF(m12,MassG)), "m12_coeff_MassG") << '\n'
            << FORMAT_ELEMENT(41, (SEMIANALYTICCOEFF(Azero,MassBp)), "Azero_coeff_MassBp") << '\n'
            << FORMAT_ELEMENT(42, (SEMIANALYTICCOEFF(m12,MassBp)), "m12_coeff_MassBp") << '\n'
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mX2m02Coeffs", SEMIANALYTICCOEFF(m02,mDx2), "m02_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2m122Coeffs", SEMIANALYTICCOEFF(m122,mDx2), "m122_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mDx2), "Azerom12_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mDx2), "Azero2_coeff_mDx2", model.get_scale());
   slha_io.set_block("mXBar2m02Coeffs", SEMIANALYTICCOEFF(m02,mDxbar2), "m02_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2m122Coeffs", SEMIANALYTICCOEFF(m122,mDxbar2), "m122_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mDxbar2), "Azerom12_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mDxbar2), "Azero2_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("TKappaAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TKappa), "Azero_coeff_TKappa", model.get_scale());
   slha_io.set_block("TKappam12Coeffs", SEMIANALYTICCOEFF(m12,TKappa), "m12_coeff_TKappa", model.get_scale());
   slha_io.set_block("TLambda12AzeroCoeffs", SEMIANALYTICCOEFF(Azero,TLambda12), "Azero_coeff_TLambda12", model.get_scale());
   slha_io.set_block("TLambda12m12Coeffs", SEMIANALYTICCOEFF(m12,TLambda12), "m12_coeff_TLambda12", model.get_scale());

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CSE6SSM_slha_io::set_extra(
   const CSE6SSM_semianalytic_slha<T>& model, const CSE6SSM_scales& scales, double MhhEFT)
{
   const CSE6SSM_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(1, (SCALES(SUSYScale)), "SUSYScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EWSBOutputParameters Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, model.get_ewsb_output_parameter(0), "m0Sqr")
            << FORMAT_ELEMENT(1, model.get_ewsb_output_parameter(1), "TanTheta")
            << FORMAT_ELEMENT(2, model.get_ewsb_output_parameter(2), "vphi")
            << FORMAT_ELEMENT(3, model.get_ewsb_output_parameter(3), "XiF")
            << FORMAT_ELEMENT(4, model.get_ewsb_output_parameter(4), "LXiF")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(2,2)/Yu(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(3,3)/Yu(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(2,2)/Yd(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(3,3)/Yd(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(2,2)/Ye(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(3,3)/Ye(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SUSYHD Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, MhhEFT, "SUSYHD Higgs pole mass")
         ;
      slha_io.set_block(block);
   }
   slha_io.set_block("TeAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYe), "Azero_coeff_TYe", model.get_scale());
   slha_io.set_block("Tem12Coeffs", SEMIANALYTICCOEFF(m12,TYe), "m12_coeff_TYe", model.get_scale());
   slha_io.set_block("TdAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYd), "Azero_coeff_TYd", model.get_scale());
   slha_io.set_block("Tdm12Coeffs", SEMIANALYTICCOEFF(m12,TYd), "m12_coeff_TYd", model.get_scale());
   slha_io.set_block("TuAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYu), "Azero_coeff_TYu", model.get_scale());
   slha_io.set_block("Tum12Coeffs", SEMIANALYTICCOEFF(m12,TYu), "m12_coeff_TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block ESIXRUNCoeffs Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(281, (SEMIANALYTICCOEFF(Azero,TKappaPr)), "Azero_coeff_TKappaPr") << '\n'
            << FORMAT_ELEMENT(282, (SEMIANALYTICCOEFF(m12,TKappaPr)), "m12_coeff_TKappaPr") << '\n'
            << FORMAT_ELEMENT(291, (SEMIANALYTICCOEFF(Azero,TSigmax)), "Azero_coeff_TSigmax") << '\n'
            << FORMAT_ELEMENT(292, (SEMIANALYTICCOEFF(m12,TSigmax)), "m12_coeff_TSigmax") << '\n'
            << FORMAT_ELEMENT(431, (SEMIANALYTICCOEFF(Azero,TSigmaL)), "Azero_coeff_TSigmaL") << '\n'
            << FORMAT_ELEMENT(432, (SEMIANALYTICCOEFF(m12,TSigmaL)), "m12_coeff_TSigmaL") << '\n'
            << FORMAT_ELEMENT(301, (SEMIANALYTICCOEFF(Azero,BMuPhi)), "Azero_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(302, (SEMIANALYTICCOEFF(m12,BMuPhi)), "m12_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(303, (SEMIANALYTICCOEFF(BMuPhi,BMuPhi)), "BMuPhi_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(304, (SEMIANALYTICCOEFF(BMuPr,BMuPhi)), "BMuPr_coeff_BMuPhi") << '\n'
            << FORMAT_ELEMENT(21, (SEMIANALYTICCOEFF(Azero,TLambdax)), "Azero_coeff_TLambdax") << '\n'
            << FORMAT_ELEMENT(22, (SEMIANALYTICCOEFF(m12,TLambdax)), "m12_coeff_TLambdax") << '\n'
            << FORMAT_ELEMENT(1011, (SEMIANALYTICCOEFF(Azero,BMuPr)), "Azero_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1012, (SEMIANALYTICCOEFF(m12,BMuPr)), "m12_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1013, (SEMIANALYTICCOEFF(BMuPhi,BMuPr)), "BMuPhi_coeff_BMuPr") << '\n'
            << FORMAT_ELEMENT(1014, (SEMIANALYTICCOEFF(BMuPr,BMuPr)), "BMuPr_coeff_BMuPr") << '\n'
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ThEAzeroCoeffs", SEMIANALYTICCOEFF(Azero,ThE), "Azero_coeff_ThE", model.get_scale());
   slha_io.set_block("ThEm12Coeffs", SEMIANALYTICCOEFF(m12,ThE), "m12_coeff_ThE", model.get_scale());
   slha_io.set_block("TgDAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TgD), "Azero_coeff_TgD", model.get_scale());
   slha_io.set_block("TgDm12Coeffs", SEMIANALYTICCOEFF(m12,TgD), "m12_coeff_TgD", model.get_scale());
   slha_io.set_block("TfuAzeroCoeffs", SEMIANALYTICCOEFF(Azero,Tfu), "Azero_coeff_Tfu", model.get_scale());
   slha_io.set_block("Tfum12Coeffs", SEMIANALYTICCOEFF(m12,Tfu), "m12_coeff_Tfu", model.get_scale());
   slha_io.set_block("TfdAzeroCoeffs", SEMIANALYTICCOEFF(Azero,Tfd), "Azero_coeff_Tfd", model.get_scale());
   slha_io.set_block("Tfdm12Coeffs", SEMIANALYTICCOEFF(m12,Tfd), "m12_coeff_Tfd", model.get_scale());
   slha_io.set_block("mq2m02Coeffs", SEMIANALYTICCOEFF(m02,mq2), "m02_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2m122Coeffs", SEMIANALYTICCOEFF(m122,mq2), "m122_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mq2), "Azerom12_coeff_mq2", model.get_scale());
   slha_io.set_block("mq2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mq2), "Azero2_coeff_mq2", model.get_scale());
   slha_io.set_block("me2m02Coeffs", SEMIANALYTICCOEFF(m02,me2), "m02_coeff_me2", model.get_scale());
   slha_io.set_block("me2m122Coeffs", SEMIANALYTICCOEFF(m122,me2), "m122_coeff_me2", model.get_scale());
   slha_io.set_block("me2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,me2), "Azerom12_coeff_me2", model.get_scale());
   slha_io.set_block("me2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,me2), "Azero2_coeff_me2", model.get_scale());
   slha_io.set_block("ml2m02Coeffs", SEMIANALYTICCOEFF(m02,ml2), "m02_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2m122Coeffs", SEMIANALYTICCOEFF(m122,ml2), "m122_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,ml2), "Azerom12_coeff_ml2", model.get_scale());
   slha_io.set_block("ml2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,ml2), "Azero2_coeff_ml2", model.get_scale());
   slha_io.set_block("mu2m02Coeffs", SEMIANALYTICCOEFF(m02,mu2), "m02_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2m122Coeffs", SEMIANALYTICCOEFF(m122,mu2), "m122_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mu2), "Azerom12_coeff_mu2", model.get_scale());
   slha_io.set_block("mu2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mu2), "Azero2_coeff_mu2", model.get_scale());
   slha_io.set_block("md2m02Coeffs", SEMIANALYTICCOEFF(m02,md2), "m02_coeff_md2", model.get_scale());
   slha_io.set_block("md2m122Coeffs", SEMIANALYTICCOEFF(m122,md2), "m122_coeff_md2", model.get_scale());
   slha_io.set_block("md2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,md2), "Azerom12_coeff_md2", model.get_scale());
   slha_io.set_block("md2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,md2), "Azero2_coeff_md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFTCoeffs Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(211, (SEMIANALYTICCOEFF(m02,mHd2)), "m02_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(212, (SEMIANALYTICCOEFF(m122,mHd2)), "m122_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(213, (SEMIANALYTICCOEFF(Azerom12,mHd2)), "Azerom12_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(214, (SEMIANALYTICCOEFF(Azero2,mHd2)), "Azero2_coeff_mHd2") << '\n'
            << FORMAT_ELEMENT(221, (SEMIANALYTICCOEFF(m02,mHu2)), "m02_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(222, (SEMIANALYTICCOEFF(m122,mHu2)), "m122_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(223, (SEMIANALYTICCOEFF(Azerom12,mHu2)), "Azerom12_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(224, (SEMIANALYTICCOEFF(Azero2,mHu2)), "Azero2_coeff_mHu2") << '\n'
            << FORMAT_ELEMENT(231, (SEMIANALYTICCOEFF(m02,ms2)), "m02_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(232, (SEMIANALYTICCOEFF(m122,ms2)), "m122_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(233, (SEMIANALYTICCOEFF(Azerom12,ms2)), "Azerom12_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(234, (SEMIANALYTICCOEFF(Azero2,ms2)), "Azero2_coeff_ms2") << '\n'
            << FORMAT_ELEMENT(241, (SEMIANALYTICCOEFF(m02,msbar2)), "m02_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(242, (SEMIANALYTICCOEFF(m122,msbar2)), "m122_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(243, (SEMIANALYTICCOEFF(Azerom12,msbar2)), "Azerom12_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(244, (SEMIANALYTICCOEFF(Azero2,msbar2)), "Azero2_coeff_msbar2") << '\n'
            << FORMAT_ELEMENT(251, (SEMIANALYTICCOEFF(m02,mphi2)), "m02_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(252, (SEMIANALYTICCOEFF(m122,mphi2)), "m122_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(253, (SEMIANALYTICCOEFF(Azerom12,mphi2)), "Azerom12_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(254, (SEMIANALYTICCOEFF(Azero2,mphi2)), "Azero2_coeff_mphi2") << '\n'
            << FORMAT_ELEMENT(261, (SEMIANALYTICCOEFF(m02,mHp2)), "m02_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(262, (SEMIANALYTICCOEFF(m122,mHp2)), "m122_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(263, (SEMIANALYTICCOEFF(Azerom12,mHp2)), "Azerom12_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(264, (SEMIANALYTICCOEFF(Azero2,mHp2)), "Azero2_coeff_mHp2") << '\n'
            << FORMAT_ELEMENT(271, (SEMIANALYTICCOEFF(m02,mHpbar2)), "m02_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(272, (SEMIANALYTICCOEFF(m122,mHpbar2)), "m122_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(273, (SEMIANALYTICCOEFF(Azerom12,mHpbar2)), "Azerom12_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(274, (SEMIANALYTICCOEFF(Azero2,mHpbar2)), "Azero2_coeff_mHpbar2") << '\n'
            << FORMAT_ELEMENT(11, (SEMIANALYTICCOEFF(Azero,MassB)), "Azero_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(12, (SEMIANALYTICCOEFF(m12,MassB)), "m12_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(21, (SEMIANALYTICCOEFF(Azero,MassWB)), "Azero_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(22, (SEMIANALYTICCOEFF(m12,MassWB)), "m12_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(31, (SEMIANALYTICCOEFF(Azero,MassG)), "Azero_coeff_MassG") << '\n'
            << FORMAT_ELEMENT(32, (SEMIANALYTICCOEFF(m12,MassG)), "m12_coeff_MassG") << '\n'
            << FORMAT_ELEMENT(41, (SEMIANALYTICCOEFF(Azero,MassBp)), "Azero_coeff_MassBp") << '\n'
            << FORMAT_ELEMENT(42, (SEMIANALYTICCOEFF(m12,MassBp)), "m12_coeff_MassBp") << '\n'
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mX2m02Coeffs", SEMIANALYTICCOEFF(m02,mDx2), "m02_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2m122Coeffs", SEMIANALYTICCOEFF(m122,mDx2), "m122_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mDx2), "Azerom12_coeff_mDx2", model.get_scale());
   slha_io.set_block("mX2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mDx2), "Azero2_coeff_mDx2", model.get_scale());
   slha_io.set_block("mXBar2m02Coeffs", SEMIANALYTICCOEFF(m02,mDxbar2), "m02_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2m122Coeffs", SEMIANALYTICCOEFF(m122,mDxbar2), "m122_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2Azerom12Coeffs", SEMIANALYTICCOEFF(Azerom12,mDxbar2), "Azerom12_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("mXBar2Azero2Coeffs", SEMIANALYTICCOEFF(Azero2,mDxbar2), "Azero2_coeff_mDxbar2", model.get_scale());
   slha_io.set_block("TKappaAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TKappa), "Azero_coeff_TKappa", model.get_scale());
   slha_io.set_block("TKappam12Coeffs", SEMIANALYTICCOEFF(m12,TKappa), "m12_coeff_TKappa", model.get_scale());
   slha_io.set_block("TLambda12AzeroCoeffs", SEMIANALYTICCOEFF(Azero,TLambda12), "Azero_coeff_TLambda12", model.get_scale());
   slha_io.set_block("TLambda12m12Coeffs", SEMIANALYTICCOEFF(m12,TLambda12), "m12_coeff_TLambda12", model.get_scale());

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CSE6SSM_slha_io::set_spectrum(const CSE6SSM<T>& model)
{
   const CSE6SSM_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CSE6SSM_slha_io::set_spectrum(const CSE6SSM_semianalytic<T>& model)
{
   const CSE6SSM_semianalytic_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void CSE6SSM_slha_io::set_spectrum(const CSE6SSM_slha<T>& model)
{
   const CSE6SSM_physical physical(model.get_physical_slha());
   const CSE6SSM_physical drbar(model.get_drbar_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_drbar_mass(drbar, model.get_scale(), write_sm_masses);
   set_drbar_mixing_matrices(drbar, model.get_scale(), write_sm_masses);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void CSE6SSM_slha_io::set_spectrum(const CSE6SSM_semianalytic_slha<T>& model)
{
   const CSE6SSM_physical physical(model.get_physical_slha());
   const CSE6SSM_physical drbar(model.get_drbar_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_drbar_mass(drbar, model.get_scale(), write_sm_masses);
   set_drbar_mixing_matrices(drbar, model.get_scale(), write_sm_masses);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#undef Pole
#undef PHYSICAL
#undef PHYSICAL_SLHA
#undef LOCALPHYSICAL
#undef MODELPARAMETER
#undef SEMIANALYTICCOEFF
#undef LowEnergyConstant
#undef SCALES

#endif
