
#ifndef DERIVS_ROOT_FINDER_H
#define DERIVS_ROOT_FINDER_H

#include <iostream>
#include <cassert>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#include "logger.hpp"
#include "error.hpp"
#include "ewsb_solver.hpp"

namespace flexiblesusy {

/**
 * @class Derivs_root_finder
 * @brief Function root finder using analytic derivatives
 *
 * The user has to provide the function (of which the root should be
 * found) of the type Function_t, and the Jacobian matrix for that
 * function of the type DFunction_t.  This function gets as arguments a
 * GSL vector of length `dimension', a pointer to the parameters (of
 * type void*) and a GSL vector where the function values must be
 * stored, or a GSL matrix where the Jacobian values must be stored.
 *
 * Example:
 * @code
 * struct Parabola {
 *    static int func(const gsl_vector* x, void*, gsl_vector* f) {
 *       const double y = gsl_vector_get(x, 0);
 *       const double z = gsl_vector_get(x, 1);
 *       gsl_vector_set(f, 0, y*(y - 5.0));
 *       gsl_vector_set(f, 1, z*(z - 1.0));
 *       return GSL_SUCCESS;
 *    }
 *
 *    static int dfunc(const gsl_vector* x, void*, gsl_matrix* J) {
 *       const double y = gsl_vector_get(x, 0);
 *       const double z = gsl_vector_get(x, 1);
 *       gsl_matrix_set(J, 0, 0, 2*y - 5.0);
 *       gsl_matrix_set(J, 0, 1, 0.0);
 *       gsl_matrix_set(J, 1, 0, 0.0);
 *       gsl_matrix_set(J, 1, 1, 2*z - 1.0);
 *       return GSL_SUCCESS;
 *    }
 *
 *    static int fdfunc(const gsl_vector* x, void*, gsl_vector* f, gsl_matrix* J) {
 *       const double y = gsl_vector_get(x, 0);
 *       const double z = gsl_vector_get(x, 1);
 *       gsl_vector_set(f, 0, y*(y - 5.0));
 *       gsl_vector_set(f, 1, z*(z - 1.0));
 *       gsl_matrix_set(J, 0, 0, 2*y - 5.0);
 *       gsl_matrix_set(J, 0, 1, 0.0);
 *       gsl_matrix_set(J, 1, 0, 0.0);
 *       gsl_matrix_set(J, 1, 1, 2*z - 1.0);
 *       return GSL_SUCCESS;
 *    }
 * };
 *
 * Derivs_root_finder<2> root_finder(Parabola::func, Parabola::dfunc,
 *                                   Parabola::fdfunc, NULL, 100, 1.0e-5);
 * const double start[2] = { 10, 10 };
 * const int status = root_finder.find_root(start);
 * @endcode
 */
template <std::size_t dimension>
class Derivs_root_finder : public EWSB_solver {
public:
   /// pointer to function to find root of
   typedef int (*Function_t)(const gsl_vector*, void*, gsl_vector*);
   /// pointer to Jacobian of function
   typedef int (*DFunction_t)(const gsl_vector*, void*, gsl_matrix*);
   /// pointer to function and Jacobian
   typedef int (*FDFunction_t)(const gsl_vector*, void*, gsl_vector*, gsl_matrix*);

   Derivs_root_finder();
   Derivs_root_finder(Function_t, DFunction_t, FDFunction_t, void*, std::size_t, 
                      double, const gsl_multiroot_fdfsolver_type* solver_type_ 
                      = gsl_multiroot_fdfsolver_hybridj);
   Derivs_root_finder(const Derivs_root_finder&);
   virtual ~Derivs_root_finder();

   double get_root(std::size_t) const;
   void set_function(Function_t f) { function = f; }
   void set_jacobian(DFunction_t j) { jacobian = j; }
   void set_combined_function(FDFunction_t c) { combined = c; }
   void set_parameters(void* m) { parameters = m; }
   void set_precision(double p) { precision = p; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   void set_solver_type(const gsl_multiroot_fdfsolver_type* t) { solver_type = t; }
   int find_root(const double[dimension]);

   // EWSB_solver interface methods
   virtual int solve(const double[dimension]);
   virtual double get_solution(unsigned);

private:
   std::size_t max_iterations; ///< maximum number of iterations
   double precision;           ///< precision goal
   gsl_vector* root;           ///< GSL vector of root
   void* parameters;           ///< pointer to parameters
   Function_t function;        ///< function to minimize
   DFunction_t jacobian;       ///< Jacobian of function
   FDFunction_t combined;      ///< function evaluating both function and Jacobian
   const gsl_multiroot_fdfsolver_type* solver_type; ///< GSL solver type

   void print_state(const gsl_multiroot_fdfsolver*, std::size_t) const;
};

/**
 * Default constructor
 */
template <std::size_t dimension>
Derivs_root_finder<dimension>::Derivs_root_finder()
   : max_iterations(100)
   , precision(1.0e-2)
   , parameters(NULL)
   , function(NULL)
   , jacobian(NULL)
   , combined(NULL)
   , solver_type(gsl_multiroot_fdfsolver_hybridj)
{
   root = gsl_vector_alloc(dimension);

   if (!root)
      throw OutOfMemoryError("GSL vector allocation failed in Derivs_root_finder()");
}

/**
 * Constructor
 *
 * @param function_ pointer to the function to minimize
 * @param jacobian_ pointer to the Jacobian of the function
 * @param parameters_ pointer to the parameters (for example the model)
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 */
template <std::size_t dimension>
Derivs_root_finder<dimension>::Derivs_root_finder(
   Function_t function_,
   DFunction_t jacobian_,
   FDFunction_t combined_,
   void* parameters_,
   std::size_t max_iterations_,
   double precision_,
   const gsl_multiroot_fdfsolver_type* solver_type_
)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , parameters(parameters_)
   , function(function_)
   , jacobian(jacobian_)
   , combined(combined_)
   , solver_type(solver_type_)
{
   root = gsl_vector_alloc(dimension);

   if (!root)
      throw OutOfMemoryError("GSL vector allocation failed in Derivs_root_finder("
                             "Function_t, DFunction_t, FDFunction_t, void*, size_t, double)");
}

template <std::size_t dimension>
Derivs_root_finder<dimension>::Derivs_root_finder(const Derivs_root_finder& other)
   : max_iterations(other.max_iterations)
   , precision(other.precision)
   , parameters(other.parameters)
   , function(other.function)
   , jacobian(other.jacobian)
   , combined(other.combined)
   , solver_type(other.solver_type)
{
   root = gsl_vector_alloc(dimension);
   gsl_vector_memcpy(root, other.root);
}

template <std::size_t dimension>
Derivs_root_finder<dimension>::~Derivs_root_finder()
{
   gsl_vector_free(root);
}

/**
 * Start the minimization
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if minimum found)
 */
template <std::size_t dimension>
int Derivs_root_finder<dimension>::find_root(const double start[dimension])
{
   assert(function && "Derivs_root_finder<dimension>::find_root: function pointer"
          " must not be zero!");
   assert(jacobian && "Derivs_root_finder<dimension>::find_root: Jacobian pointer"
          " must not be zero!");
   assert(combined && "Derivs_root_finder<dimension>::find_root: Combined pointer"
          " must not be zero!");

   int status;
   std::size_t iter = 0;
   gsl_multiroot_function_fdf f = {function, jacobian, combined,
                                   dimension, parameters};

   gsl_multiroot_fdfsolver* solver
      = gsl_multiroot_fdfsolver_alloc(solver_type, dimension);

   if (!solver) {
      throw OutOfMemoryError(std::string("Cannot allocate gsl_multiroot_fdfsolver ") +
                             gsl_multiroot_fdfsolver_name(solver));
   }

#ifndef ENABLE_DEBUG
   gsl_set_error_handler_off();
#endif

   for (std::size_t i = 0; i < dimension; ++i)
      gsl_vector_set(root, i, start[i]);

   gsl_multiroot_fdfsolver_set(solver, &f, root);

#ifdef ENABLE_VERBOSE
   print_state(solver, iter);
#endif

   do {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate(solver);

#ifdef ENABLE_VERBOSE
      print_state(solver, iter);
#endif

      if (status)   // check if solver is stuck
         break;

      status = gsl_multiroot_test_residual(solver->f, precision);
   } while (status == GSL_CONTINUE && iter < max_iterations);

#ifdef ENABLE_VERBOSE
   printf("\tRoot_finder status = %s\n", gsl_strerror(status));
#endif

   gsl_vector_memcpy(root, solver->x);
   gsl_multiroot_fdfsolver_free(solver);

   return status;
}

/**
 * Print state of the root finder
 *
 * @param solver solver
 * @param iteration iteration number
 */
template <std::size_t dimension>
void Derivs_root_finder<dimension>::print_state(const gsl_multiroot_fdfsolver* solver,
                                         std::size_t iteration) const
{
   std::cout << "\tIteration " << iteration << ": x =";
   for (std::size_t i = 0; i < dimension; ++i)
      std::cout << ' ' << gsl_vector_get(solver->x, i);
   std::cout << ", f(x) =";
   for (std::size_t i = 0; i < dimension; ++i)
      std::cout << ' ' << gsl_vector_get(solver->f, i);
   std::cout << '\n';
}

template <std::size_t dimension>
double Derivs_root_finder<dimension>::get_root(std::size_t i) const
{
   assert(i < dimension && "Derivs_root_finder<>::get_root: index out"
          " of bounds");
   return gsl_vector_get(root, i);
}

template <std::size_t dimension>
int Derivs_root_finder<dimension>::solve(const double start[dimension])
{
   return (find_root(start) == GSL_SUCCESS ?
           EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

template <std::size_t dimension>
double Derivs_root_finder<dimension>::get_solution(unsigned i)
{
   return get_root(i);
}

} // namespace flexiblesusy

#endif
