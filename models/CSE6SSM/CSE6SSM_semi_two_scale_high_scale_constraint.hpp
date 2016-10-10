// ====================================================================
// Test implementation of high scale constraint to be used in
// semianalytic version of the two-scale algorithm
// ====================================================================

#ifndef CSE6SSM_SEMI_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define CSE6SSM_SEMI_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "CSE6SSM_semi_high_scale_constraint.hpp"
#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CSE6SSM_semianalytic;

class Two_scale;

template<>
class CSE6SSM_semianalytic_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CSE6SSM_semianalytic_high_scale_constraint();
   CSE6SSM_semianalytic_high_scale_constraint(CSE6SSM_semianalytic<Two_scale>*);
   virtual ~CSE6SSM_semianalytic_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CSE6SSM_semianalytic_input_parameters<Two_scale>& get_input_parameters() const;
   CSE6SSM_semianalytic<Two_scale>* get_model() const;
   void initialize();
   void set_initial_scale_guess(double); ///< set initial scale guess before iteration
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CSE6SSM_semianalytic<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
