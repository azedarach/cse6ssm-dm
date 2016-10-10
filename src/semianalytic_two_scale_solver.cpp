// ====================================================================
// Test implementation of semianalytic solver routines, using the
// two-scale algorithm to solve the inner iteration
// ====================================================================

#include "semianalytic_two_scale_solver.hpp"
#include "two_scale_constraint.hpp"
#include "two_scale_convergence_tester.hpp"
#include "two_scale_initial_guesser.hpp"
#include "two_scale_matching.hpp"
#include "two_scale_model.hpp"
#include "two_scale_running_precision.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "functors.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cassert>

/**
 * @file semianalytic_two_scale_solver.cpp
 * @brief contains the implementation of the RGFlow_semianalytic<Two_scale>
 */

namespace flexiblesusy {

RGFlow_semianalytic<Two_scale>::RGFlow_semianalytic()
   : models()
   , iteration(0)
   , inner_convergence_tester(NULL)
   , outer_convergence_tester(NULL)
   , initial_guesser(NULL)
   , running_precision_calculator(NULL)
   , running_precision(1.0e-3)
   , inner_solver()
   , model_at_this_scale(NULL)
{
}

RGFlow_semianalytic<Two_scale>::~RGFlow_semianalytic()
{
   delete_models();
}

/**
 * @brief Solves the boundary value problem.
 *
 * At first the initial_guess() is called.  Afterwards, the semianalytic
 * solution is applied: first, the function
 * iteratively runs the tower up and down and imposes the boundary
 * conditions.  The iteration stops if either the maximum number of
 * iterations is reached or the precision goal is achieved (defined by
 * the convergence_tester for the inner two-scale iteration). Then
 * the next, semianalytic step is applied. This is iterated until
 * either the maximum number of iterations is reached or the precision goal
 * is achieved (defined by the convergence_tester for the outer iteration).
 */
void RGFlow_semianalytic<Two_scale>::solve()
{
   check_setup();

   const unsigned max_iterations = get_max_iterations();
   if (models.empty() || max_iterations == 0)
      return;

   initial_guess();

   iteration = 0;
   bool accuracy_reached = false;
   const size_t number_of_models = models.size();
   while (iteration < max_iterations && !accuracy_reached) {
      update_running_precision();
      clear_problems();

      // solve inner iteration
      inner_convergence_tester->restart();

      inner_solver.reset();
      inner_solver.set_convergence_tester(inner_convergence_tester);
      inner_solver.set_running_precision(running_precision_calculator);
      for (size_t m = 0; m < number_of_models; ++m) {
         const TModel* model = models[m];
         inner_solver.add_model(model->model, model->matching_condition,
                                model->inner_upwards_constraints,
                                model->inner_downwards_constraints);
      }
      inner_solver.solve();

      run_semianalytic_up();
      run_semianalytic_down();
      accuracy_reached = accuracy_goal_reached();
      ++iteration;
   }

   apply_lowest_constraint();

   if (!accuracy_reached)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

/**
 * Sanity checks the models and boundary conditions
 */
void RGFlow_semianalytic<Two_scale>::check_setup() const
{
   for (size_t m = 0, num_models = models.size(); m < num_models; ++m) {
      const TModel* model = models[m];
      if (!model->model) {
         std::stringstream message;
         message << "RGFlow_semianalytic<Two_scale>::Error: model pointer ["
                 << m << "] is NULL";
         throw SetupError(message.str());
      }

      // check whether last model has a non-zero matching condition
      if (m + 1 == num_models) {
         if (model->matching_condition)
            WARNING("the matching condition of the " << model->model->name()
                    << " is non-zero but will not be used");
      } else {
         if (model->matching_condition == NULL) {
            std::stringstream message;
            message << "RGFlow_semianalytic<Two_scale>::Error: matching condition "
                    << "of the " << model->model->name() << " to the "
                    << models[m + 1]->model->name() << " is NULL";
            throw SetupError(message.str());
         }
      }
   }

   if (!inner_convergence_tester) {
      throw SetupError("RGFlow_semianalytic<Two_scale>::Error: inner convergence tester must "
                       "not be NULL");
   }
   if (!outer_convergence_tester) {
      throw SetupError("RGFlow_semianalytic<Two_scale>::Error: outer convergence tester must "
                       "not be NULL");
   }
}

void RGFlow_semianalytic<Two_scale>::clear_problems()
{
   VERBOSE_MSG("> clearing problems ...");

   const size_t number_of_models = models.size();
   for (size_t m = 0; m < number_of_models; ++m) {
      TModel* model = models[m];
      model->model->clear_problems();
   }
}

void RGFlow_semianalytic<Two_scale>::delete_models()
{
   for_each(models.begin(), models.end(), Delete_object());
}

/**
 * Does the initial guess by calling the guess() method of the initial
 * guesser (if given)
 */
void RGFlow_semianalytic<Two_scale>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

/**
 * Run the model tower to the highest scale.  Thereby all upwards semianalytic
 * constraints are imposed (by calling the apply() functions) and the
 * matching conditions are applied (by calling
 * match_low_to_high_scale_model() functions).
 */
void RGFlow_semianalytic<Two_scale>::run_semianalytic_up()
{
   VERBOSE_MSG(">> running tower up (iteration " << iteration << ") ...");
   const size_t number_of_models = models.size();
   for (size_t m = 0; m < number_of_models; ++m) {
      TModel* model = models[m];
      model->model->set_precision(get_precision());
      VERBOSE_MSG(">> \tselecting model " << model->model->name());
      // apply all constraints
      const size_t n_upwards_constraints = model->semianalytic_upwards_constraints.size();
      for (size_t c = 0; c < n_upwards_constraints; ++c) {
         Constraint<Two_scale>* constraint = model->semianalytic_upwards_constraints[c];
         const double scale = constraint->get_scale();
         VERBOSE_MSG(">> \t\tselecting constraint " << c << " at scale " << scale);
         VERBOSE_MSG(">> \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         VERBOSE_MSG(">> \t\t\tapplying constraint");
         constraint->apply();
      }
      // apply matching condition if this is not the last model
      if (m != number_of_models - 1) {
         VERBOSE_MSG(">> \tmatching to model " << models[m + 1]->model->name());
         Matching<Two_scale>* mc = model->matching_condition;
         const double scale = mc->get_scale();
         VERBOSE_MSG(">> \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         VERBOSE_MSG(">> \t\t\tapplying matching condition");
         mc->match_low_to_high_scale_model();
      }
   }
   VERBOSE_MSG(">> running up finished");
}

/**
 * Run the model tower to the lowest scale.  Thereby all downwards semianalytic
 * constraints are imposed (by calling the apply() functions) and the
 * matching conditions are applied (by calling
 * match_high_to_low_scale_model() functions).
 */
void RGFlow_semianalytic<Two_scale>::run_semianalytic_down()
{
   assert(models.size() > 0 && "model size must not be zero");
   VERBOSE_MSG("<< running tower down ...");
   const size_t number_of_models = models.size();
   for (size_t m = number_of_models; m-- > 0;) {
      TModel* model = models[m];
      VERBOSE_MSG("<< \tselecting model " << model->model->name());
      // apply all semianalytic constraints
      // If m is the last model, do not apply the highest constraint
      // because it was already applied when we ran up
      const size_t c_begin = (m + 1 == number_of_models ? 1 : 0);
      const size_t c_end = model->semianalytic_downwards_constraints.size();
      for (size_t c = c_begin; c < c_end; ++c) {
         Constraint<Two_scale>* constraint = model->semianalytic_downwards_constraints[c];
         const double scale = constraint->get_scale();
         VERBOSE_MSG("<< \t\trunning model to scale " << scale);
         model->model->run_to(scale);
         // If m is the lowest energy model, do not apply the lowest
         // constraint, because it will be applied when we run up next
         // time.
         if (m != 0 || c + 1 != c_end) {
            VERBOSE_MSG("<< \t\t\tapply constraint");
            constraint->apply();
         }
      }
      // apply matching condition if this is not the first model
      if (m > 0) {
         Matching<Two_scale>* mc = models[m - 1]->matching_condition;
         VERBOSE_MSG("<< \tmatching to model " << models[m - 1]->model->name());
         const double scale = mc->get_scale();
         VERBOSE_MSG("<< \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         VERBOSE_MSG("<< \t\t\tapplying matching condition");
         mc->match_high_to_low_scale_model();
      }
   }
   VERBOSE_MSG("<< running down finished");
}

/**
 * Impose the constraint at the lowest energy scale (by calling the
 * apply() function).
 */
void RGFlow_semianalytic<Two_scale>::apply_lowest_constraint()
{
   if (models.empty())
      return;

   TModel* model = models[0];
   model_at_this_scale = model->model;

   // determine whether lowest scale occurs in inner or
   // semianalytic iteration
   double inner_lowest_scale, semianalytic_lowest_scale;

   if (model->inner_downwards_constraints.empty())
      inner_lowest_scale = 0.;
   else
      inner_lowest_scale = model->inner_downwards_constraints.back()->get_scale();

   if (model->semianalytic_downwards_constraints.empty())
      semianalytic_lowest_scale = 0.;
   else
      semianalytic_lowest_scale = model->semianalytic_downwards_constraints.back()->get_scale();

   Constraint<Two_scale>* constraint = (inner_lowest_scale < semianalytic_lowest_scale ? 
                                        model->inner_downwards_constraints.back() :
                                        model->semianalytic_downwards_constraints.back());
   const double scale = constraint->get_scale();
   VERBOSE_MSG("| selecting constraint 0 at scale " << scale);
   VERBOSE_MSG("| \trunning model " << model->model->name() << " to scale " << scale);
   model->model->run_to(scale);
   VERBOSE_MSG("| \tapplying constraint");
   constraint->apply();
}

/**
 * Returns the precision of the RG running
 *
 * @return RG running precision
 */
double RGFlow_semianalytic<Two_scale>::get_precision()
{
   return running_precision;
}

/**
 * Recalculates the precision of the RG running using the user defined
 * Two_scale_running_precision class.
 */
void RGFlow_semianalytic<Two_scale>::update_running_precision()
{
   if (running_precision_calculator)
      running_precision = running_precision_calculator->get_precision(iteration);
}

/**
 * Add a model and the corresponding model constraints. Note that the
 * order of the model registration is important: Models that are added
 * later are assumed to be valid at a higher scale. The same is true
 * for the constraints: they are assumed to be ordered from low to
 * high energies.
 *
 * @param model model
 * @param inner_constraints vector of model constraints for inner iteration
 * @param semianalytic_constraints vector of model constraints for outer iteration
 */
void RGFlow_semianalytic<Two_scale>::add_model(Two_scale_model* model,
                                               const std::vector<Constraint<Two_scale>*>& inner_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_constraints)
{
   add_model(model, NULL, inner_constraints, semianalytic_constraints);
}

/**
 * Add a model and the corresponding model constraints.  With this
 * function the user can use different constraints for the up and down
 * running of the model tower.
 *
 * @param model model
 * @param inner_upwards_constraints model constraints for running up in inner iteration
 * @param inner_downwards_constraints model constraints for running down in inner iteration
 * @param semianalytic_upwards_constraints model constraints for running up in outer iteration
 * @param semianalytic_downwards_constraints model constraints for running down in outer iteration
 */
void RGFlow_semianalytic<Two_scale>::add_model(Two_scale_model* model,
                                               const std::vector<Constraint<Two_scale>*>& inner_upwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& inner_downwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_upwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_downwards_constraints)
{
   add_model(model, NULL, inner_upwards_constraints, inner_downwards_constraints,
             semianalytic_upwards_constraints, semianalytic_downwards_constraints);
}

/**
 * Add a model, the corresponding model constraints and the matching
 * condition to the next model.
 * @param model model
 * @param mc matching condition to the next higher model
 * @param inner_constraints vector of model constraints for inner iteration
 * @param semianalytic_constraints vector of model constraints for outer iteration
 */
void RGFlow_semianalytic<Two_scale>::add_model(Two_scale_model* model,
                                               Matching<Two_scale>* mc,
                                               const std::vector<Constraint<Two_scale>*>& inner_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_constraints)
{
   // create vectors of downward constraints
   std::vector<Constraint<Two_scale>*> inner_downwards_constraints;
   std::reverse_copy(inner_constraints.begin(), inner_constraints.end(),
                     std::back_inserter(inner_downwards_constraints));
   std::vector<Constraint<Two_scale>*> semianalytic_downwards_constraints;
   std::reverse_copy(semianalytic_constraints.begin(), semianalytic_constraints.end(),
                     std::back_inserter(semianalytic_downwards_constraints));

   add_model(model, mc, inner_constraints, inner_downwards_constraints,
             semianalytic_constraints, semianalytic_downwards_constraints);
}

/**
 * Add a model, the corresponding model constraints and the matching
 * condition to the next higher model.  With this function the user
 * can use different constraints for the up and down running of the
 * model tower.
 *
 * @param model model
 * @param mc matching condition to the next higher model
 * @param inner_upwards_constraints model constraints for running up in inner iteration
 * @param inner_downwards_constraints model constraints for running down in inner iteration
 * @param semianalytic_upwards_constraints model constraints for running up in outer iteration
 * @param semianalytic_downwards_constraints model constraints for running down in outer iteration
 */
void RGFlow_semianalytic<Two_scale>::add_model(Two_scale_model* model,
                                               Matching<Two_scale>* mc,
                                               const std::vector<Constraint<Two_scale>*>& inner_upwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& inner_downwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_upwards_constraints,
                                               const std::vector<Constraint<Two_scale>*>& semianalytic_downwards_constraints)
{
   TModel* new_model = new TModel(model, inner_upwards_constraints, inner_downwards_constraints,
                                  semianalytic_upwards_constraints, semianalytic_downwards_constraints, mc);

   for (Constraint_container::iterator it = new_model->inner_upwards_constraints.begin(),
           end = new_model->inner_upwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->inner_downwards_constraints.begin(),
           end = new_model->inner_downwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->semianalytic_upwards_constraints.begin(),
           end = new_model->semianalytic_upwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->semianalytic_downwards_constraints.begin(),
           end = new_model->semianalytic_downwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   if (!models.empty())
      models.back()->matching_condition->set_models(models.back()->model,
                                                    model);

   models.push_back(new_model);
}

/**
 * Returns the value returned by the accuracy_goal_reached() method of
 * the convergence tester.
 */
bool RGFlow_semianalytic<Two_scale>::accuracy_goal_reached() const
{
   return outer_convergence_tester->accuracy_goal_reached();
}

/**
 * Set the convergence testers to be used during the iteration.
 *
 * @param inner_tester the convergence tester to be used in the inner iteration
 * @param outer_tester the convergence tester to be used in the outer iteration
 */
void RGFlow_semianalytic<Two_scale>::set_convergence_testers(Convergence_tester<Two_scale>* inner_tester, Convergence_tester<Two_scale>* outer_tester)
{
   inner_convergence_tester = inner_tester;
   outer_convergence_tester = outer_tester;
}

/**
 * Set the initial guesser.
 *
 * @param ig initial guesser
 */
void RGFlow_semianalytic<Two_scale>::set_initial_guesser(Initial_guesser<Two_scale>* ig)
{
   initial_guesser = ig;
}

/**
 * Set RG running precision calculator.
 *
 * @param rp running precision calculator
 */
void RGFlow_semianalytic<Two_scale>::set_running_precision(Two_scale_running_precision* rp)
{
   running_precision_calculator = rp;
}

/**
 * Returns the number of performed iterations
 * @return number of performed iterations
 */
unsigned int RGFlow_semianalytic<Two_scale>::number_of_iterations_done() const
{
   return iteration;
}

/**
 * Returns the maximum number of iterations set in the convergence
 * tester.
 */
unsigned int RGFlow_semianalytic<Two_scale>::get_max_iterations() const
{
   return outer_convergence_tester->max_iterations();
}

/**
 * @brief resets the solver to the initial condition
 *
 * The pointers to the models, matching conditions, convergence
 * tester, initial guesser, and running precision calculator are set
 * to zero.  The running precision is set to the default value 0.001.
 */
void RGFlow_semianalytic<Two_scale>::reset()
{
   delete_models();
   models.clear();
   inner_solver.reset();

   iteration = 0;
   inner_convergence_tester = NULL;
   outer_convergence_tester = NULL;
   initial_guesser = NULL;
   running_precision_calculator = NULL;
   running_precision = 1.0e-3;
   model_at_this_scale = NULL;
}

/**
 * Run the model tower to the given scale.
 *
 * @param scale scale to run to
 */
void RGFlow_semianalytic<Two_scale>::run_to(double scale)
{
   // find model which is defined at `scale'
   model_at_this_scale = NULL;
   const size_t number_of_models = models.size();

   for (size_t m = 0; m < number_of_models; ++m) {
      const TModel* model = models[m];
      double highest_scale, lowest_scale;

      if (!model) {
         std::ostringstream msg;
         msg << "RGFlow_semianalytic<Two_scale>::run_to: pointer to model " << m
             << " is zero";
         throw SetupError(msg.str());
      }

      if (m != number_of_models - 1) {
         // if this is not the last model, the matching condition is
         // the highest scale
         const Matching<Two_scale>* mc = model->matching_condition;
         if (!mc) {
            std::ostringstream msg;
            msg << "RGFlow_semianalytic<Two_scale>::run_to: pointer to matching condition"
                  " of model " << m << " is zero";
            throw SetupError(msg.str());
         }
         highest_scale = mc->get_scale();
      } else {
         // otherwise the last constraint is at the highest scale
         double inner_highest_scale, semianalytic_highest_scale;

         if (model->inner_upwards_constraints.empty())
            inner_highest_scale = std::numeric_limits<double>::max();
         else
            inner_highest_scale = model->inner_upwards_constraints.back()->get_scale();

         if (model->semianalytic_upwards_constraints.empty())
            semianalytic_highest_scale = std::numeric_limits<double>::max();
         else
            semianalytic_highest_scale = model->semianalytic_upwards_constraints.back()->get_scale();

         highest_scale = (inner_highest_scale > semianalytic_highest_scale ?
                          inner_highest_scale : semianalytic_highest_scale);
      }

      if (m > 0) {
         // if this is not the first model, the previous matching
         // condition is the lowest scale
         lowest_scale = models[m-1]->matching_condition->get_scale();
      } else {
         // otherwise the first constraint is at the lowest scale
         double inner_lowest_scale, semianalytic_lowest_scale;

         if (model->inner_upwards_constraints.empty())
            inner_lowest_scale = 0.;
         else
            inner_lowest_scale = model->inner_upwards_constraints[0]->get_scale();

         if (model->semianalytic_upwards_constraints.empty())
            semianalytic_lowest_scale = 0.;
         else
            semianalytic_lowest_scale = model->semianalytic_upwards_constraints[0]->get_scale();

         lowest_scale = (inner_lowest_scale < semianalytic_lowest_scale ?
                         inner_lowest_scale : semianalytic_lowest_scale);
      }

      if (lowest_scale <= scale && scale <= highest_scale) {
         model_at_this_scale = model->model;
         break;
      }
   }

   if (model_at_this_scale)
      model_at_this_scale->run_to(scale);
   else {
      std::ostringstream msg;
      msg << "No model at scale " << scale << " found!";
      throw SetupError(msg.str());
   }
}

/**
 * Returns the pointer to the model at the current scale.
 */
Two_scale_model* RGFlow_semianalytic<Two_scale>::get_model() const
{
   return model_at_this_scale;
}

} // namespace flexiblesusy
