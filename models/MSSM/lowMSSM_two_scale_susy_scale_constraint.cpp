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

// File generated at Wed 15 Jul 2015 13:12:13

#include "lowMSSM_two_scale_susy_scale_constraint.hpp"
#include "lowMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"

#include <cassert>
#include <cmath>

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
#define MODELCLASSNAME lowMSSM<Two_scale>

lowMSSM_susy_scale_constraint<Two_scale>::lowMSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

lowMSSM_susy_scale_constraint<Two_scale>::lowMSSM_susy_scale_constraint(
   lowMSSM<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
{
   initialize();
}

lowMSSM_susy_scale_constraint<Two_scale>::~lowMSSM_susy_scale_constraint()
{
}

void lowMSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: lowMSSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto TYdInput = INPUTPARAMETER(TYdInput);
   const auto TYeInput = INPUTPARAMETER(TYeInput);
   const auto TYuInput = INPUTPARAMETER(TYuInput);
   const auto BMuInput = INPUTPARAMETER(BMuInput);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto MassBInput = INPUTPARAMETER(MassBInput);
   const auto MassWBInput = INPUTPARAMETER(MassWBInput);
   const auto MassGInput = INPUTPARAMETER(MassGInput);

   MODEL->set_Mu(Re(MuInput));
   MODEL->set_TYd((TYdInput).real());
   MODEL->set_TYe((TYeInput).real());
   MODEL->set_TYu((TYuInput).real());
   MODEL->set_BMu(Re(BMuInput));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_MassB(Re(MassBInput));
   MODEL->set_MassWB(Re(MassWBInput));
   MODEL->set_MassG(Re(MassGInput));


   // the parameters, which are fixed by the EWSB eqs., will now be
   // defined at this scale (at the EWSB loop level defined in the
   // model)
   model->solve_ewsb();
}

double lowMSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double lowMSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const lowMSSM_input_parameters<Two_scale>& lowMSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   assert(model && "Error: lowMSSM_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

lowMSSM<Two_scale>* lowMSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void lowMSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<lowMSSM<Two_scale>*>(model_);
}

void lowMSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void lowMSSM_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "lowMSSM_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 1000;

   scale = initial_scale_guess;
}

void lowMSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "lowMSSM_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const auto ZU = MODELPARAMETER(ZU);
   const auto MSu = MODELPARAMETER(MSu);

   scale = Sqrt(Power(MSu(0),Sqr(Abs(ZU(0,2))) + Sqr(Abs(ZU(0,5))))*Power(MSu(1
      ),Sqr(Abs(ZU(1,2))) + Sqr(Abs(ZU(1,5))))*Power(MSu(2),Sqr(Abs(ZU(2,2))) +
      Sqr(Abs(ZU(2,5))))*Power(MSu(3),Sqr(Abs(ZU(3,2))) + Sqr(Abs(ZU(3,5))))*Power
      (MSu(4),Sqr(Abs(ZU(4,2))) + Sqr(Abs(ZU(4,5))))*Power(MSu(5),Sqr(Abs(ZU(5,2))
      ) + Sqr(Abs(ZU(5,5)))));


}

} // namespace flexiblesusy
