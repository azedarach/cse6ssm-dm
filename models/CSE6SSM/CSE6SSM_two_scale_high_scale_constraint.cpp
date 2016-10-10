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

// File generated at Wed 3 Jun 2015 23:47:49

#include "CSE6SSM_two_scale_high_scale_constraint.hpp"
#include "CSE6SSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "numerics2.hpp"

#include <cassert>
#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CSE6SSM<Two_scale>

CSE6SSM_high_scale_constraint<Two_scale>::CSE6SSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

CSE6SSM_high_scale_constraint<Two_scale>::CSE6SSM_high_scale_constraint(
   CSE6SSM<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
{
   initialize();
}

CSE6SSM_high_scale_constraint<Two_scale>::~CSE6SSM_high_scale_constraint()
{
}

void CSE6SSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: CSE6SSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CSE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(3.54491);
   }
   if (std::fabs(model->get_g2()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CSE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(3.54491);
   }
   if (std::fabs(model->get_g3()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CSE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(3.54491);
   }

   update_scale();

   const auto MuPrInput = INPUTPARAMETER(MuPrInput);
   const auto BMuPrInput = INPUTPARAMETER(BMuPrInput);
   const auto MuPhiInput = INPUTPARAMETER(MuPhiInput);
   const auto BMuPhiInput = INPUTPARAMETER(BMuPhiInput);
   const auto SigmaLInput = INPUTPARAMETER(SigmaLInput);
   const auto KappaPrInput = INPUTPARAMETER(KappaPrInput);
   const auto SigmaxInput = INPUTPARAMETER(SigmaxInput);
   const auto gDInput = INPUTPARAMETER(gDInput);
   const auto hEInput = INPUTPARAMETER(hEInput);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto Lambda12Input = INPUTPARAMETER(Lambda12Input);
   const auto fuInput = INPUTPARAMETER(fuInput);
   const auto fdInput = INPUTPARAMETER(fdInput);
   const auto Azero = INPUTPARAMETER(Azero);
   const auto m0 = INPUTPARAMETER(m0);
   const auto m12 = INPUTPARAMETER(m12);
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   MODEL->set_g1p(Re(g1));
   MODEL->set_MuPr(Re(MuPrInput));
   MODEL->set_BMuPr(Re(BMuPrInput));
   MODEL->set_MuPhi(Re(MuPhiInput));
   MODEL->set_BMuPhi(Re(BMuPhiInput));
   MODEL->set_SigmaL(Re(SigmaLInput));
   MODEL->set_KappaPr(Re(KappaPrInput));
   MODEL->set_Sigmax(Re(SigmaxInput));
   MODEL->set_gD((gDInput).real());
   MODEL->set_hE((hEInput).real());
   MODEL->set_Kappa((KappaInput).real());
   MODEL->set_Lambda12((Lambda12Input).real());
   MODEL->set_fu((fuInput).real());
   MODEL->set_fd((fdInput).real());
   MODEL->set_TYe((Azero*Ye).real());
   MODEL->set_TYd((Azero*Yd).real());
   MODEL->set_TYu((Azero*Yu).real());
   MODEL->set_TKappaPr(Re(Azero*KappaPrInput));
   MODEL->set_TSigmax(Re(Azero*SigmaxInput));
   MODEL->set_ThE((Azero*hEInput).real());
   MODEL->set_TSigmaL(Re(Azero*SigmaLInput));
   MODEL->set_TgD((Azero*gDInput).real());
   MODEL->set_Tfu((Azero*fuInput).real());
   MODEL->set_Tfd((Azero*fdInput).real());
   MODEL->set_TKappa((Azero*KappaInput).real());
   MODEL->set_TLambda12((Azero*Lambda12Input).real());
   MODEL->set_TLambdax(Re(Azero*Lambdax));
   MODEL->set_mHd2(Re(Sqr(m0)));
   MODEL->set_mHu2(Re(Sqr(m0)));
   MODEL->set_ms2(Re(Sqr(m0)));
   MODEL->set_msbar2(Re(Sqr(m0)));
   MODEL->set_mphi2(Re(Sqr(m0)));
   MODEL->set_mHp2(Re(Sqr(m0)));
   MODEL->set_mHpbar2(Re(Sqr(m0)));
   MODEL->set_mH1I2((Sqr(m0)*UNITMATRIX(2)).real());
   MODEL->set_mH2I2((Sqr(m0)*UNITMATRIX(2)).real());
   MODEL->set_mSI2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_mq2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_ml2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_md2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_mu2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_me2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_mDx2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_mDxbar2((Sqr(m0)*UNITMATRIX(3)).real());
   MODEL->set_MassB(Re(m12));
   MODEL->set_MassWB(Re(m12));
   MODEL->set_MassG(Re(m12));
   MODEL->set_MassBp(Re(m12));

   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto g1p = MODELPARAMETER(g1p);
      const auto Yd = MODELPARAMETER(Yd);
      const auto hE = MODELPARAMETER(hE);
      const auto Ye = MODELPARAMETER(Ye);
      const auto SigmaL = MODELPARAMETER(SigmaL);
      const auto KappaPr = MODELPARAMETER(KappaPr);
      const auto Sigmax = MODELPARAMETER(Sigmax);
      const auto gD = MODELPARAMETER(gD);
      const auto Kappa = MODELPARAMETER(Kappa);
      const auto Lambda12 = MODELPARAMETER(Lambda12);
      const auto Lambdax = MODELPARAMETER(Lambdax);
      const auto fu = MODELPARAMETER(fu);
      const auto fd = MODELPARAMETER(fd);
      const auto Yu = MODELPARAMETER(Yu);

      if (MaxAbsValue(g1) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g1", MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g1");
      if (MaxAbsValue(g2) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g2", MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g2");
      if (MaxAbsValue(g3) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g3", MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g3");
      if (MaxAbsValue(g1p) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g1p", MaxAbsValue(g1p), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g1p");
      if (MaxAbsValue(Yd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yd");
      if (MaxAbsValue(hE) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("hE", MaxAbsValue(hE), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("hE");
      if (MaxAbsValue(Ye) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Ye");
      if (MaxAbsValue(SigmaL) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("SigmaL", MaxAbsValue(SigmaL), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("SigmaL");
      if (MaxAbsValue(KappaPr) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("KappaPr", MaxAbsValue(KappaPr), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("KappaPr");
      if (MaxAbsValue(Sigmax) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Sigmax", MaxAbsValue(Sigmax), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Sigmax");
      if (MaxAbsValue(gD) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("gD", MaxAbsValue(gD), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("gD");
      if (MaxAbsValue(Kappa) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Kappa", MaxAbsValue(Kappa), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Kappa");
      if (MaxAbsValue(Lambda12) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Lambda12", MaxAbsValue(Lambda12), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Lambda12");
      if (MaxAbsValue(Lambdax) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Lambdax", MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Lambdax");
      if (MaxAbsValue(fu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("fu", MaxAbsValue(fu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("fu");
      if (MaxAbsValue(fd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("fd", MaxAbsValue(fd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("fd");
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yu");

   }
}

double CSE6SSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double CSE6SSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const CSE6SSM_input_parameters<Two_scale>& CSE6SSM_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

CSE6SSM<Two_scale>* CSE6SSM_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void CSE6SSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<CSE6SSM<Two_scale>*>(model_);
}

void CSE6SSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void CSE6SSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void CSE6SSM_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "CSE6SSM_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 2.e16;

   scale = initial_scale_guess;
}

void CSE6SSM_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "CSE6SSM_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const CSE6SSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*Exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("CSE6SSM_high_scale_constraint<Two_scale>: Overflow error"
            " during calculation of high scale: " << strerror(errno) << '\n'
            << "   current scale = " << currentScale << '\n'
            << "   new scale = " << scale << '\n'
            << "   resetting scale to " << get_initial_scale_guess());
#endif
      scale = get_initial_scale_guess();
      errno = 0;
   }


}

} // namespace flexiblesusy
