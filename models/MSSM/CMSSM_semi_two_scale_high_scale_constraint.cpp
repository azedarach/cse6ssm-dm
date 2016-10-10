
#include "CMSSM_semi_two_scale_high_scale_constraint.hpp"
#include "CMSSM_semi_two_scale_model.hpp"
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
#define MODELCLASSNAME CMSSM_semianalytic<Two_scale>

CMSSM_semianalytic_high_scale_constraint<Two_scale>::CMSSM_semianalytic_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

CMSSM_semianalytic_high_scale_constraint<Two_scale>::CMSSM_semianalytic_high_scale_constraint(CMSSM_semianalytic<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
{
   initialize();
}

CMSSM_semianalytic_high_scale_constraint<Two_scale>::~CMSSM_semianalytic_high_scale_constraint()
{
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: CMSSM_semianalytic_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(3.54491);
   }
   if (std::fabs(model->get_g2()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(3.54491);
   }
   if (std::fabs(model->get_g3()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(3.54491);
   }

   update_scale();

   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto MuInput_at_MS = INPUTPARAMETER(MuInput_at_MS);

   if (!MuInput_at_MS) {
      MODEL->set_Mu(Re(MuInput));
   }

   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto Yd = MODELPARAMETER(Yd);
      const auto Ye = MODELPARAMETER(Ye);
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
      if (MaxAbsValue(Yd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yd");
      if (MaxAbsValue(Ye) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Ye");
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yu");
   }
}

double CMSSM_semianalytic_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double CMSSM_semianalytic_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const CMSSM_semianalytic_input_parameters<Two_scale>& CMSSM_semianalytic_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

CMSSM_semianalytic<Two_scale>* CMSSM_semianalytic_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<CMSSM_semianalytic<Two_scale>*>(model_);
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::set_initial_scale_guess(double s)
{
   initial_scale_guess = s;
   scale = initial_scale_guess;
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "CMSSM_semianalytic_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 2.e16;

   scale = initial_scale_guess;
}

void CMSSM_semianalytic_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "CMSSM_semianalytic_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const MSSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*Exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("CMSSM_semianalytic_high_scale_constraint<Two_scale>: Overflow error"
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
