
#ifndef CMSSM_SEMI_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define CMSSM_SEMI_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "CMSSM_semi_susy_scale_constraint.hpp"
#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "CMSSM_semi_two_scale_high_scale_constraint.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CMSSM_semianalytic;

class Two_scale;

template<>
class CMSSM_semianalytic_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CMSSM_semianalytic_susy_scale_constraint();
   CMSSM_semianalytic_susy_scale_constraint(CMSSM_semianalytic<Two_scale>*);
   virtual ~CMSSM_semianalytic_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CMSSM_semianalytic_input_parameters<Two_scale>& get_input_parameters() const;
   CMSSM_semianalytic<Two_scale>* get_model() const;
   void initialize();

   CMSSM_semianalytic_high_scale_constraint<Two_scale>* get_input_scale_constraint() const;
   void set_input_scale_constraint(CMSSM_semianalytic_high_scale_constraint<Two_scale>*);

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CMSSM_semianalytic<Two_scale>* model;
   CMSSM_semianalytic_high_scale_constraint<Two_scale>* high_constraint;
};

} // namespace flexiblesusy

#endif
