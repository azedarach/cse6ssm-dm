// ====================================================================
// Specialisation of a constraint handler for setting up the solver
// when using the two_scale algorithm.
// ====================================================================

#ifndef CSE6SSM_TWO_SCALE_CONSTRAINT_HANDLER_H
#define CSE6SSM_TWO_SCALE_CONSTRAINT_HANDLER_H

#include "CSE6SSM_constraint_handler.hpp" 
#include "CSE6SSM_two_scale_high_scale_constraint.hpp"
#include "CSE6SSM_two_scale_susy_scale_constraint.hpp"
#include "CSE6SSM_two_scale_low_scale_constraint.hpp"
#include "CSE6SSM_two_scale_initial_guesser.hpp"

#include "two_scale_solver.hpp"

namespace flexiblesusy {

class Two_scale;
class Two_scale_model;

template<>
class CSE6SSM_constraint_handler<Two_scale> {
public:
   CSE6SSM_constraint_handler();
   virtual ~CSE6SSM_constraint_handler();

   double get_highest_scale() const;
   double get_susy_scale() const;
   double get_lowest_scale() const;

   void initialize_constraints(Two_scale_model*, const QedQcd&);
   CSE6SSM_initial_guesser<Two_scale> get_initial_guesser(CSE6SSM<Two_scale>*,
                                                          const QedQcd&) const;
   void add_constraints_to_solver(Two_scale_model*, RGFlow<Two_scale>&);

   void set_sm_parameters(const QedQcd& oneset);
   const QedQcd& get_sm_parameters() const;

private:
   CSE6SSM_high_scale_constraint<Two_scale> high_scale_constraint;
   CSE6SSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   CSE6SSM_low_scale_constraint<Two_scale> low_scale_constraint;
};

} // namespace flexiblesusy

#endif
