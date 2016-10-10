
#include "CMSSM_semi_two_scale_constraint_handler.hpp"

namespace flexiblesusy {

CMSSM_semianalytic_constraint_handler<Two_scale>::CMSSM_semianalytic_constraint_handler()
   : high_scale_constraint()
   , susy_scale_constraint()
   , low_scale_constraint()
{
}

CMSSM_semianalytic_constraint_handler<Two_scale>::~CMSSM_semianalytic_constraint_handler()
{
}

double CMSSM_semianalytic_constraint_handler<Two_scale>::get_highest_scale() const
{
   return high_scale_constraint.get_scale();
}

double CMSSM_semianalytic_constraint_handler<Two_scale>::get_susy_scale() const
{
   return susy_scale_constraint.get_scale();
}

double CMSSM_semianalytic_constraint_handler<Two_scale>::get_lowest_scale() const
{
   return low_scale_constraint.get_scale();
}

void CMSSM_semianalytic_constraint_handler<Two_scale>::initialize_constraints(Two_scale_model* model, const QedQcd& oneset)
{
   high_scale_constraint.clear();
   susy_scale_constraint.clear();
   low_scale_constraint.clear();

   high_scale_constraint.set_model(model);
   susy_scale_constraint.set_model(model);
   low_scale_constraint.set_model(model);

   low_scale_constraint.set_sm_parameters(oneset);

   susy_scale_constraint.set_input_scale_constraint(&high_scale_constraint);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();
}

CMSSM_semianalytic_initial_guesser<Two_scale> CMSSM_semianalytic_constraint_handler<Two_scale>::get_initial_guesser(CMSSM_semianalytic<Two_scale>* model, const QedQcd& oneset)
{
   CMSSM_semianalytic_initial_guesser<Two_scale> initial_guesser(model, oneset,
                                                                 low_scale_constraint,
                                                                 susy_scale_constraint,
                                                                 high_scale_constraint,
                                                                 &high_scale_constraint);

   return initial_guesser;
}

void CMSSM_semianalytic_constraint_handler<Two_scale>::add_constraints_to_solver(Two_scale_model* model, RGFlow_semianalytic<Two_scale>& solver)
{
   std::vector<Constraint<Two_scale>*> inner_upward_constraints(2);
   inner_upward_constraints[0] = &low_scale_constraint;
   inner_upward_constraints[1] = &high_scale_constraint;

   std::vector<Constraint<Two_scale>*> inner_downward_constraints(2);
   inner_downward_constraints[0] = &high_scale_constraint;
   inner_downward_constraints[1] = &low_scale_constraint;

   std::vector<Constraint<Two_scale>*> outer_upward_constraints(1);
   outer_upward_constraints[0] = &susy_scale_constraint;

   std::vector<Constraint<Two_scale>*> outer_downward_constraints(2);
   outer_downward_constraints[0] = &susy_scale_constraint;
   outer_downward_constraints[1] = &low_scale_constraint;

   solver.add_model(model, inner_upward_constraints, inner_downward_constraints,
                    outer_upward_constraints, outer_downward_constraints);
}

void CMSSM_semianalytic_constraint_handler<Two_scale>::set_sm_parameters(const QedQcd& oneset)
{
   low_scale_constraint.set_sm_parameters(oneset);
}

const QedQcd& CMSSM_semianalytic_constraint_handler<Two_scale>::get_sm_parameters() const
{
   return low_scale_constraint.get_sm_parameters();
}

} // namespace flexiblesusy
