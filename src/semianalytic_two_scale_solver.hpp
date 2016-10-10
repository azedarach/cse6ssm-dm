// ====================================================================
// Test implementation of semianalytic solver routines, using the
// two-scale algorithm to solve the inner iteration
// ====================================================================

#ifndef SEMIANALYTIC_TWO_SCALE_SOLVER_H
#define SEMIANALYTIC_TWO_SCALE_SOLVER_H

#include "rg_flow_semianalytic.hpp"
#include "two_scale_solver.hpp"

#include <vector>
#include <string>
#include <sstream>

/**
 * @file semianalytic_two_scale_solver.hpp
 * @brief contains the definition of the RGFlow_semianalytic<Two_scale> class
 */

namespace flexiblesusy {

template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
class Two_scale;
class Two_scale_model;
class Two_scale_running_precision;

/**
 * @class RGFlow_semianalytic<Two_scale>
 * @brief Boundary condition solver (semianalytic two-scale algorithm)
 *
 * This boundary condition solver uses a semianalytic algorithm
 * to solve the boundary value problem. A two-scale iteration is
 * used to determine an initial set of parameters, followed by
 * a second step in which these parameters are used to to derive
 * the remaining parameters using a semi-analytic approach provided
 * by the user. The whole process is iterated until convergence is
 * reached.
 */

template<>
class RGFlow_semianalytic<Two_scale> {
public:
   RGFlow_semianalytic();
   ~RGFlow_semianalytic();

   /// add models and constraints
   void add_model(Two_scale_model*,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&);
   /// add models and constraints
   void add_model(Two_scale_model*,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&);
   /// add model, constraints and matching condition
   void add_model(Two_scale_model*,
                  Matching<Two_scale>* m = NULL,
                  const std::vector<Constraint<Two_scale>*>&
                  inner_constraints = std::vector<Constraint<Two_scale>*>(),
                  const std::vector<Constraint<Two_scale>*>&
                  semianalytic_constraints = std::vector<Constraint<Two_scale>*>());
   /// add models, constraints and matching condition
   void add_model(Two_scale_model*,
                  Matching<Two_scale>*,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&);
   /// get model at current scale
   Two_scale_model* get_model() const;
   /// get number of used iterations
   unsigned int number_of_iterations_done() const;
   /// clear all internal data
   void reset();
   /// pick valid model and run it to the given scale
   void run_to(double);
   /// set convergence testers
   void set_convergence_testers(Convergence_tester<Two_scale>*,
                                Convergence_tester<Two_scale>*);
   /// set running precision calculator
   void set_running_precision(Two_scale_running_precision*);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser<Two_scale>*);
   /// solves the boundary value problem
   void solve();

private:
   typedef std::vector<Constraint<Two_scale>*> Constraint_container;

   /**
    * @class TModel
    * @brief contains model, constraints and matching condition
    *
    * This class lumps together the model, its constraints and the
    * matching condition to the next higher model.
    */
   struct TModel {
      Two_scale_model* model;                          ///< the model
      Constraint_container inner_upwards_constraints;  ///< upwards two-scale constraints
      Constraint_container inner_downwards_constraints;  ///< downwards two-scale constraints
      Constraint_container semianalytic_upwards_constraints; ///< upwards semianalytic constraints
      Constraint_container semianalytic_downwards_constraints; ///< downwards semianalytic constraints
      Matching<Two_scale>* matching_condition;         ///< matching condition

      TModel(Two_scale_model* m,
             const Constraint_container& iuc,
             const Constraint_container& idc,
             const Constraint_container& suc,
             const Constraint_container& sdc,
             Matching<Two_scale>* mc)
         : model(m)
         , inner_upwards_constraints(iuc)
         , inner_downwards_constraints(idc)
         , semianalytic_upwards_constraints(suc)
         , semianalytic_downwards_constraints(sdc)
         , matching_condition(mc)
         {}
   };

   std::vector<TModel*> models; ///< tower of models (from low to high scale)
   unsigned int iteration; ///< iteration number (starting at 0)
   Convergence_tester<Two_scale>* inner_convergence_tester; ///< the convergence tester for the two-scale iteration
   Convergence_tester<Two_scale>* outer_convergence_tester; ///< the convergence tester for the overall iteration
   Initial_guesser<Two_scale>* initial_guesser; ///< does initial guess
   Two_scale_running_precision* running_precision_calculator; ///< RG running precision calculator
   double running_precision; ///< RG running precision used
   RGFlow<Two_scale> inner_solver; ///< solves inner two-scale iteration
   Two_scale_model* model_at_this_scale; ///< model at current scale

   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const; ///< check the setup
   void clear_problems(); ///< clear model problems
   void delete_models(); ///< delete all models
   unsigned int get_max_iterations() const; /// < returns max. number of iterations
   void initial_guess(); ///< initial guess
   void run_semianalytic_up(); ///< run through upwards semianalytic constraints
   void run_semianalytic_down(); ///< run through downwards semianalytic constraints
   void apply_lowest_constraint(); ///< apply lowest constraint
   double get_precision(); ///< returns running precision
   void update_running_precision(); ///< update the RG running precision
};

} // namespace flexiblesusy

#endif
