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

#ifndef CMSSM_SLHA_IO_H
#define CMSSM_SLHA_IO_H

#include "CMSSM_semi_two_scale_model_slha.hpp"
#include "CMSSM_two_scale_model_slha.hpp"
#include "MSSM_info.hpp"
#include "MSSM_physical.hpp"
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
struct CMSSM_input_parameters;

template <class T>
struct CMSSM_semianalytic_input_parameters;

class Two_scale;

class Spectrum_generator_settings;

struct CMSSM_scales {
   CMSSM_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class CMSSM_slha_io {
public:
   CMSSM_slha_io();
   ~CMSSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(CMSSM_input_parameters<Two_scale>&) const;
   void fill(CMSSM_semianalytic_input_parameters<Two_scale>&) const;
   void fill(MSSM_mass_eigenstates&) const;
   template <class T> void fill(CMSSM_slha<T>&) const;
   template <class T> void fill(CMSSM_semianalytic_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const CMSSM_input_parameters<Two_scale>&);
   void set_extpar(const CMSSM_semianalytic_input_parameters<Two_scale>&);
   template <class T> void set_extra(const CMSSM_slha<T>&, const CMSSM_scales&);
   template <class T> void set_extra(const CMSSM_semianalytic_slha<T>&, const CMSSM_scales&);
   void set_minpar(const CMSSM_input_parameters<Two_scale>&);
   void set_minpar(const CMSSM_semianalytic_input_parameters<Two_scale>&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const CMSSM_slha<T>&);
   template <class T> void set_spectrum(const CMSSM<T>&);
   template <class T> void set_spectrum(const CMSSM_semianalytic_slha<T>&);
   template <class T> void set_spectrum(const CMSSM_semianalytic<T>&);
   void set_spinfo(const Problems<MSSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(CMSSM_input_parameters<Two_scale>&, int, double);
   static void fill_extpar_tuple(CMSSM_input_parameters<Two_scale>&, int, double);
   static void fill_minpar_tuple(CMSSM_semianalytic_input_parameters<Two_scale>&, int, double);
   static void fill_extpar_tuple(CMSSM_semianalytic_input_parameters<Two_scale>&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CMSSM_slha<T>&, const QedQcd&, const CMSSM_scales&);
   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CMSSM_semianalytic_slha<T>&, const QedQcd&, const CMSSM_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CMSSM_slha<T>&, const QedQcd&);
   template <class T>
   static SLHAea::Coll fill_slhaea(const CMSSM_semianalytic_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CMSSM_slha<T>&, const QedQcd&, const CMSSM_scales&);
   template <class T>
   static SLHAea::Coll fill_slhaea(const CMSSM_semianalytic_slha<T>&, const QedQcd&, const CMSSM_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 14;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_drbar_mass(const MSSM_physical&, double, bool);
   void set_drbar_mixing_matrices(const MSSM_physical&, double, bool);
   void set_mass(const MSSM_physical&, bool);
   void set_mixing_matrices(const MSSM_physical&, bool);
   template <class T> void set_model_parameters(const CMSSM_slha<T>&);
   template <class T> void set_model_parameters(const CMSSM_semianalytic_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(MSSM_mass_eigenstates&) const;
   void fill_drbar(MSSM_physical&) const;
   void fill_physical(MSSM_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CMSSM_slha_io::fill(CMSSM_slha<T>& model) const
{
   fill(static_cast<MSSM_mass_eigenstates&>(model));
   fill_drbar(model.get_drbar_slha());
   fill_physical(model.get_physical_slha());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CMSSM_slha_io::fill(CMSSM_semianalytic_slha<T>& model) const
{
   fill(static_cast<MSSM_mass_eigenstates&>(model));
   fill_drbar(model.get_drbar_slha());
   fill_physical(model.get_physical_slha());
}

template <class T>
void CMSSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CMSSM_slha<T>& model,
   const QedQcd& qedqcd, const CMSSM_scales& scales)
{
   CMSSM_slha_io slha_io;
   const CMSSM_input_parameters<Two_scale>& input = model.get_input();
   const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems
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
void CMSSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CMSSM_semianalytic_slha<T>& model,
   const QedQcd& qedqcd, const CMSSM_scales& scales)
{
   CMSSM_slha_io slha_io;
   const CMSSM_semianalytic_input_parameters<Two_scale>& input = model.get_input();
   const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems
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
SLHAea::Coll CMSSM_slha_io::fill_slhaea(
   const CMSSM_slha<T>& model, const QedQcd& qedqcd)
{
   CMSSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CMSSM_slha_io::fill_slhaea(
   const CMSSM_semianalytic_slha<T>& model, const QedQcd& qedqcd)
{
   CMSSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CMSSM_slha_io::fill_slhaea(
   const CMSSM_slha<T>& model, const QedQcd& qedqcd,
   const CMSSM_scales& scales)
{
   SLHAea::Coll slhaea;
   CMSSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

template <class T>
SLHAea::Coll CMSSM_slha_io::fill_slhaea(
   const CMSSM_semianalytic_slha<T>& model, const QedQcd& qedqcd,
   const CMSSM_scales& scales)
{
   SLHAea::Coll slhaea;
   CMSSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CMSSM_slha_io::set_model_parameters(const CMSSM_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
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
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
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
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
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
void CMSSM_slha_io::set_model_parameters(const CMSSM_semianalytic_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
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
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
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
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
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
void CMSSM_slha_io::set_extra(
   const CMSSM_slha<T>& model, const CMSSM_scales& scales)
{
   const MSSM_physical physical(model.get_physical_slha());

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
      block << "Block ALPHA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_NUMBER((ArcSin(Pole(ZH(1,1)))), "ArcSin(Pole(ZH(1,1)))")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(vu)/MODELPARAMETER(vd)), "vu/vd")
            << FORMAT_ELEMENT(3, (Sqrt(Sqr(MODELPARAMETER(vd)) + Sqr(MODELPARAMETER(vu)))), "Sqrt(Sqr(vd) + Sqr(vu))")
            << FORMAT_ELEMENT(4, (Sqr(MODELPARAMETER(MAh)(1))), "Sqr(MAh(1))")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(0,0)/Yu(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(2,2)/Yu(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(0,0)/Yd(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(2,2)/Yd(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(0,0)/Ye(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(2,2)/Ye(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(31, (Sqrt(MODELPARAMETER(ml2)(0,0))), "Sqrt(ml2(0,0))")
            << FORMAT_ELEMENT(32, (Sqrt(MODELPARAMETER(ml2)(1,1))), "Sqrt(ml2(1,1))")
            << FORMAT_ELEMENT(33, (Sqrt(MODELPARAMETER(ml2)(2,2))), "Sqrt(ml2(2,2))")
            << FORMAT_ELEMENT(34, (Sqrt(MODELPARAMETER(me2)(0,0))), "Sqrt(me2(0,0))")
            << FORMAT_ELEMENT(35, (Sqrt(MODELPARAMETER(me2)(1,1))), "Sqrt(me2(1,1))")
            << FORMAT_ELEMENT(36, (Sqrt(MODELPARAMETER(me2)(2,2))), "Sqrt(me2(2,2))")
            << FORMAT_ELEMENT(41, (Sqrt(MODELPARAMETER(mq2)(0,0))), "Sqrt(mq2(0,0))")
            << FORMAT_ELEMENT(42, (Sqrt(MODELPARAMETER(mq2)(1,1))), "Sqrt(mq2(1,1))")
            << FORMAT_ELEMENT(43, (Sqrt(MODELPARAMETER(mq2)(2,2))), "Sqrt(mq2(2,2))")
            << FORMAT_ELEMENT(44, (Sqrt(MODELPARAMETER(mu2)(0,0))), "Sqrt(mu2(0,0))")
            << FORMAT_ELEMENT(45, (Sqrt(MODELPARAMETER(mu2)(1,1))), "Sqrt(mu2(1,1))")
            << FORMAT_ELEMENT(46, (Sqrt(MODELPARAMETER(mu2)(2,2))), "Sqrt(mu2(2,2))")
            << FORMAT_ELEMENT(47, (Sqrt(MODELPARAMETER(md2)(0,0))), "Sqrt(md2(0,0))")
            << FORMAT_ELEMENT(48, (Sqrt(MODELPARAMETER(md2)(1,1))), "Sqrt(md2(1,1))")
            << FORMAT_ELEMENT(49, (Sqrt(MODELPARAMETER(md2)(2,2))), "Sqrt(md2(2,2))")
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
void CMSSM_slha_io::set_extra(
   const CMSSM_semianalytic_slha<T>& model, const CMSSM_scales& scales)
{
   const MSSM_physical physical(model.get_physical_slha());

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
            << FORMAT_ELEMENT(1, model.get_ewsb_output_parameter(1), "BMu0")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block ALPHA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_NUMBER((ArcSin(Pole(ZH(1,1)))), "ArcSin(Pole(ZH(1,1)))")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(vu)/MODELPARAMETER(vd)), "vu/vd")
            << FORMAT_ELEMENT(3, (Sqrt(Sqr(MODELPARAMETER(vd)) + Sqr(MODELPARAMETER(vu)))), "Sqrt(Sqr(vd) + Sqr(vu))")
            << FORMAT_ELEMENT(4, (Sqr(MODELPARAMETER(MAh)(1))), "Sqr(MAh(1))")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(0,0)/Yu(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(2,2)/Yu(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(0,0)/Yd(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(2,2)/Yd(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(0,0)/Ye(0,0)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(2,2)/Ye(2,2)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(31, (Sqrt(MODELPARAMETER(ml2)(0,0))), "Sqrt(ml2(0,0))")
            << FORMAT_ELEMENT(32, (Sqrt(MODELPARAMETER(ml2)(1,1))), "Sqrt(ml2(1,1))")
            << FORMAT_ELEMENT(33, (Sqrt(MODELPARAMETER(ml2)(2,2))), "Sqrt(ml2(2,2))")
            << FORMAT_ELEMENT(34, (Sqrt(MODELPARAMETER(me2)(0,0))), "Sqrt(me2(0,0))")
            << FORMAT_ELEMENT(35, (Sqrt(MODELPARAMETER(me2)(1,1))), "Sqrt(me2(1,1))")
            << FORMAT_ELEMENT(36, (Sqrt(MODELPARAMETER(me2)(2,2))), "Sqrt(me2(2,2))")
            << FORMAT_ELEMENT(41, (Sqrt(MODELPARAMETER(mq2)(0,0))), "Sqrt(mq2(0,0))")
            << FORMAT_ELEMENT(42, (Sqrt(MODELPARAMETER(mq2)(1,1))), "Sqrt(mq2(1,1))")
            << FORMAT_ELEMENT(43, (Sqrt(MODELPARAMETER(mq2)(2,2))), "Sqrt(mq2(2,2))")
            << FORMAT_ELEMENT(44, (Sqrt(MODELPARAMETER(mu2)(0,0))), "Sqrt(mu2(0,0))")
            << FORMAT_ELEMENT(45, (Sqrt(MODELPARAMETER(mu2)(1,1))), "Sqrt(mu2(1,1))")
            << FORMAT_ELEMENT(46, (Sqrt(MODELPARAMETER(mu2)(2,2))), "Sqrt(mu2(2,2))")
            << FORMAT_ELEMENT(47, (Sqrt(MODELPARAMETER(md2)(0,0))), "Sqrt(md2(0,0))")
            << FORMAT_ELEMENT(48, (Sqrt(MODELPARAMETER(md2)(1,1))), "Sqrt(md2(1,1))")
            << FORMAT_ELEMENT(49, (Sqrt(MODELPARAMETER(md2)(2,2))), "Sqrt(md2(2,2))")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("TeAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYe), "Azero_coeff_TYe", model.get_scale());
   slha_io.set_block("Tem12Coeffs", SEMIANALYTICCOEFF(m12,TYe), "m12_coeff_TYe", model.get_scale());
   slha_io.set_block("TdAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYd), "Azero_coeff_TYd", model.get_scale());
   slha_io.set_block("Tdm12Coeffs", SEMIANALYTICCOEFF(Azero,TYd), "m12_coeff_TYd", model.get_scale());
   slha_io.set_block("TuAzeroCoeffs", SEMIANALYTICCOEFF(Azero,TYu), "Azero_coeff_TYu", model.get_scale());
   slha_io.set_block("Tum12Coeffs", SEMIANALYTICCOEFF(Azero,TYu), "m12_coeff_TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIXCoeffs Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1011, (SEMIANALYTICCOEFF(Azero,BMu)), "Azero_coeff_BMu") << '\n'
            << FORMAT_ELEMENT(1012, (SEMIANALYTICCOEFF(m12,BMu)), "m12_coeff_BMu") << '\n'
            << FORMAT_ELEMENT(1013, (SEMIANALYTICCOEFF(BMu0,BMu)), "BMu0_coeff_BMu") << '\n'
         ;
      slha_io.set_block(block);
   }
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
            << FORMAT_ELEMENT(11, (SEMIANALYTICCOEFF(Azero,MassB)), "Azero_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(12, (SEMIANALYTICCOEFF(m12,MassB)), "m12_coeff_MassB") << '\n'
            << FORMAT_ELEMENT(21, (SEMIANALYTICCOEFF(Azero,MassWB)), "Azero_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(22, (SEMIANALYTICCOEFF(m12,MassWB)), "m12_coeff_MassWB") << '\n'
            << FORMAT_ELEMENT(31, (SEMIANALYTICCOEFF(Azero,MassG)), "Azero_coeff_MassG") << '\n'
            << FORMAT_ELEMENT(32, (SEMIANALYTICCOEFF(m12,MassG)), "m12_coeff_MassG") << '\n'
         ;
      slha_io.set_block(block);
   }
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CMSSM_slha_io::set_spectrum(const CMSSM<T>& model)
{
   const CMSSM_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CMSSM_slha_io::set_spectrum(const CMSSM_semianalytic<T>& model)
{
   const CMSSM_semianalytic_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void CMSSM_slha_io::set_spectrum(const CMSSM_slha<T>& model)
{
   const MSSM_physical physical(model.get_physical_slha());
   const MSSM_physical drbar(model.get_drbar_slha());
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
void CMSSM_slha_io::set_spectrum(const CMSSM_semianalytic_slha<T>& model)
{
   const MSSM_physical physical(model.get_physical_slha());
   const MSSM_physical drbar(model.get_drbar_slha());
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
#undef LowEnergyConstant
#undef SCALES

#endif
